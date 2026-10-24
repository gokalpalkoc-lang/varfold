"""Live folding: validation, per-client rate limiting and the Tailscale worker client."""

from __future__ import annotations

import ipaddress
import re
import time
from collections import defaultdict, deque

import httpx
import numpy as np

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class LiveError(Exception):
    def __init__(self, status: int, detail: str, retry_after: int | None = None):
        super().__init__(detail)
        self.status, self.detail, self.retry_after = status, detail, retry_after


def clean_sequence(seq: str, max_length: int) -> str:
    seq = re.sub(r"\s+", "", seq).upper()
    if not seq:
        raise LiveError(422, "empty sequence")
    if seq.startswith(">"):
        raise LiveError(422, "paste the bare sequence without a FASTA header")
    bad = sorted(set(seq) - AMINO_ACIDS)
    if bad:
        raise LiveError(422, f"unsupported residue letters: {''.join(bad)}")
    if len(seq) > max_length:
        raise LiveError(422, f"sequence has {len(seq)} residues; the live demo accepts up to {max_length}")
    return seq


def point_mutations(reference: str, seq: str, offset: int) -> list[str] | None:
    """Substitutions relative to the reference construct, or None if not comparable."""
    if len(seq) != len(reference):
        return None
    return [f"{a}{offset + i}{b}" for i, (a, b) in enumerate(zip(reference, seq)) if a != b]


def client_key(host: str | None) -> str:
    """Rate-limit identity: the IPv4 address, or the /64 for IPv6 (one subscriber's block)."""
    if not host:
        return "unknown"
    try:
        addr = ipaddress.ip_address(host)
    except ValueError:
        return host
    if addr.version == 6:
        if addr.ipv4_mapped is not None:
            return str(addr.ipv4_mapped)
        return str(ipaddress.ip_network(f"{addr}/64", strict=False))
    return str(addr)


class RateLimiter:
    """Sliding-window limit per client key (single-process; fine for one uvicorn worker)."""

    def __init__(self, limit: int, window_s: float = 3600.0):
        self.limit, self.window = limit, window_s
        self._hits: dict[str, deque] = defaultdict(deque)

    def check(self, key: str, now: float | None = None) -> None:
        now = time.monotonic() if now is None else now
        hits = self._hits[key]
        while hits and now - hits[0] > self.window:
            hits.popleft()
        if len(hits) >= self.limit:
            retry = int(self.window - (now - hits[0])) + 1
            raise LiveError(429, f"live-fold limit reached ({self.limit}/hour); "
                                 f"precomputed variants stay instant", retry_after=retry)
        hits.append(now)
        if len(self._hits) > 10_000:  # bound memory: forget clients whose last hit has expired
            for k in [k for k, v in self._hits.items() if not v or now - v[-1] > self.window]:
                del self._hits[k]

    def refund(self, key: str) -> None:
        """Give back the most recent hit (the fold never ran: worker busy/offline)."""
        if hits := self._hits.get(key):
            hits.pop()


class WorkerClient:
    def __init__(self, base_url: str, token: str, timeout_s: float, transport: httpx.AsyncBaseTransport | None = None):
        self._client = httpx.AsyncClient(base_url=base_url.rstrip("/"), transport=transport,
                                         headers={"Authorization": f"Bearer {token}"},
                                         timeout=httpx.Timeout(timeout_s, connect=5.0))
        self._status: tuple[float, dict] | None = None

    async def close(self) -> None:
        await self._client.aclose()

    async def status(self, max_age_s: float = 10.0) -> dict:
        if self._status and time.monotonic() - self._status[0] < max_age_s:
            return self._status[1]
        try:
            resp = await self._client.get("/health", timeout=4.0)
            resp.raise_for_status()
            h = resp.json()
            status = {"available": True, "busy": h["busy"], "queue": h["queue"], "max_queue": h["max_queue"],
                      "max_length": h["max_length"], "device": h["device"], "fingerprint": h["fingerprint"]}
        except (httpx.HTTPError, KeyError, ValueError):
            status = {"available": False}
        self._status = (time.monotonic(), status)
        return status

    async def fold(self, sequence: str) -> dict:
        try:
            resp = await self._client.post("/fold", json={"sequence": sequence})
        except httpx.ConnectTimeout as exc:  # a TimeoutException, but it means "unreachable"
            self._status = None
            raise LiveError(503, "the live engine is offline right now; precomputed variants still work") from exc
        except httpx.TimeoutException as exc:
            raise LiveError(504, "the live engine took too long; try a shorter sequence") from exc
        except httpx.HTTPError as exc:
            self._status = None
            raise LiveError(503, "the live engine is offline right now; precomputed variants still work") from exc
        if resp.status_code == 503:
            raise LiveError(503, "the live engine is busy; try again in a minute", retry_after=30)
        if resp.status_code == 422:
            raise LiveError(422, resp.json().get("detail", "invalid sequence"))
        if resp.status_code != 200:
            raise LiveError(502, "the live engine returned an error")
        return resp.json()


# --------------------------------------------------------------------------- geometry

def _ca(pdb: str) -> np.ndarray:
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in pdb.splitlines()
                     if l.startswith("ATOM") and l[12:16].strip() == "CA"])


def superpose_onto(pdb: str, reference_pdb: str) -> str:
    """Rigidly move ``pdb`` onto ``reference_pdb`` (equal length, Cα Kabsch)."""
    mobile, target = _ca(pdb), _ca(reference_pdb)
    if len(mobile) != len(target) or len(mobile) < 3:
        return pdb
    mc, tc = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = tc - mc @ r.T
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            x = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ r.T + t
            line = f"{line[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{line[54:]}"
        out.append(line)
    return "\n".join(out) + "\n"
