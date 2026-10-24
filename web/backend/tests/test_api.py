"""API tests against a bundle produced by the real exporter (``varfold export``).

Requires the engine package (installed in the repo venv) to build the fixture.
"""

import gzip
import json
import sys
from pathlib import Path

import httpx
import pytest
from fastapi.testclient import TestClient

ENGINE_TESTS = Path(__file__).resolve().parents[3] / "engine" / "tests"
sys.path.insert(0, str(ENGINE_TESTS))

pytest.importorskip("varfold")
from synthetic import build_synthetic_run  # noqa: E402

from app.main import create_app  # noqa: E402
from app.settings import Settings  # noqa: E402

TOKEN = "t" * 32


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    from varfold.analyze import run_analysis
    from varfold.export import export_web_db

    cfg = build_synthetic_run(tmp_path_factory.mktemp("run"))
    run_analysis(cfg)
    return Path(export_web_db(cfg)), cfg


class FakeWorker:
    """Stands in for the GPU worker behind Tailscale."""

    def __init__(self, bundle_path: Path):
        self.calls = 0
        self.online = True
        self.busy = False
        import sqlite3
        con = sqlite3.connect(bundle_path)
        self.pdb = gzip.decompress(con.execute("SELECT pdb_gz FROM variants WHERE variant_id='WT'").fetchone()[0]).decode()
        self.fingerprint = json.loads(dict(con.execute("SELECT key, value FROM meta"))["model"])["fingerprint"]

    def handler(self, request: httpx.Request) -> httpx.Response:
        if not self.online:
            raise httpx.ConnectError("down", request=request)
        if request.headers.get("authorization") != f"Bearer {TOKEN}":
            return httpx.Response(401)
        if request.url.path == "/health":
            return httpx.Response(200, json={"status": "ok", "busy": False, "queue": 0, "max_queue": 3,
                                             "max_length": 400, "device": "fake", "fingerprint": self.fingerprint})
        if self.busy:
            return httpx.Response(503, json={"detail": "busy"})
        self.calls += 1
        seq = json.loads(request.content)["sequence"]
        return httpx.Response(200, json={"pdb": self.pdb, "plddt": [80.0] * len(seq), "mean_plddt": 80.0,
                                         "ptm": 0.8, "runtime_s": 1.0, "cached": False, "device": "fake",
                                         "fingerprint": self.fingerprint})


@pytest.fixture
def client_and_worker(bundle):
    path, _ = bundle
    worker = FakeWorker(path)
    settings = Settings(db_path=path, static_dir=None, worker_url="http://worker:8100", worker_token=TOKEN,
                        live_per_hour=2)
    app = create_app(settings, worker_transport=httpx.MockTransport(worker.handler))
    with TestClient(app) as client:
        yield client, worker


def test_meta_and_listing(client_and_worker, bundle):
    client, _ = client_and_worker
    meta = client.get("/api/meta").json()
    assert meta["construct"]["start"] == 101
    assert len(meta["wild_type"]["plddt"]) == 80
    assert meta["live"]["enabled"] is True

    res = client.get("/api/variants", params={"label": "pathogenic,benign", "sort": "delta_local_plddt"}).json()
    assert res["total"] == 26
    values = [i["delta_local_plddt"] for i in res["items"]]
    assert values == sorted(values)
    assert res["items"][0]["label"] == "pathogenic"

    breaking = client.get("/api/variants", params={"breaking": True, "limit": 500}).json()
    assert breaking["total"] >= 14 and all(i["breaking_point"] == 1 for i in breaking["items"])


def test_listing_rejects_bad_input(client_and_worker):
    client, _ = client_and_worker
    assert client.get("/api/variants", params={"sort": "pdb_gz; DROP TABLE variants"}).status_code == 422
    assert client.get("/api/variants", params={"label": "nope"}).status_code == 422
    assert client.get("/api/variants", params={"limit": 10_000}).status_code == 422
    # LIKE wildcards in the search box are literal, not patterns.
    assert client.get("/api/variants", params={"q": "%"}).json()["total"] == 0


def test_variant_detail_and_structure(client_and_worker):
    client, _ = client_and_worker
    first = client.get("/api/variants", params={"limit": 1}).json()["items"][0]
    v = client.get(f"/api/variants/{first['variant_id']}").json()
    assert len(v["plddt"]) == 80 and v["sequence"]
    assert client.get("/api/variants/NOPE1A").status_code == 404

    r = client.get(f"/api/variants/{first['variant_id']}/structure.pdb")
    assert r.status_code == 200 and "ATOM" in r.text
    etag = r.headers["etag"]
    assert client.get(f"/api/variants/{first['variant_id']}/structure.pdb",
                      headers={"If-None-Match": etag}).status_code == 304
    assert len(client.get("/api/positions").json()) == 80


def test_live_fold_dataset_hit_skips_worker_and_limit(client_and_worker, bundle):
    client, worker = client_and_worker
    _, cfg = bundle
    seq = json.loads(cfg.paths.dataset_meta.read_text())["construct"]["sequence"]
    for _ in range(5):  # well beyond the per-hour limit
        r = client.post("/api/live/fold", json={"sequence": seq.lower()})
        assert r.status_code == 200 and r.json() == {"source": "dataset", "variant_id": "WT",
                                                     "sequence": seq, "mutations": []}
    assert worker.calls == 0


def test_live_fold_worker_path_and_rate_limit(client_and_worker, bundle):
    client, worker = client_and_worker
    _, cfg = bundle
    seq = json.loads(cfg.paths.dataset_meta.read_text())["construct"]["sequence"]
    novel = ("W" if seq[0] != "W" else "Y") + ("W" if seq[1] != "W" else "Y") + seq[2:]
    r = client.post("/api/live/fold", json={"sequence": novel})
    body = r.json()
    assert r.status_code == 200 and body["source"] == "live" and body["comparable"] is True
    assert len(body["mutations"]) == 2 and body["mutations"][0].endswith(novel[0])
    assert client.post("/api/live/fold", json={"sequence": "ACDEFGHIK"}).status_code == 200
    limited = client.post("/api/live/fold", json={"sequence": "ACDEFGHIKL"})
    assert limited.status_code == 429 and "Retry-After" in limited.headers
    assert worker.calls == 2


def test_live_fold_errors(client_and_worker):
    client, worker = client_and_worker
    assert client.post("/api/live/fold", json={"sequence": "ACDXZ"}).status_code == 422
    assert client.post("/api/live/fold", json={"sequence": "A" * 401}).status_code == 422
    worker.busy = True
    assert client.post("/api/live/fold", json={"sequence": "MKTAYIAK"}).status_code == 503
    worker.busy, worker.online = False, False
    assert client.post("/api/live/fold", json={"sequence": "MKTAYIAKQ"}).status_code == 503
    assert client.get("/api/live/status").json()["available"] is False


def test_live_status(client_and_worker):
    client, _ = client_and_worker
    s = client.get("/api/live/status").json()
    assert s["available"] is True and s["configured"] is True and s["max_length"] == 400


def test_spa_fallback_and_traversal(bundle, tmp_path):
    path, _ = bundle
    static = tmp_path / "dist"
    (static / "assets").mkdir(parents=True)
    (static / "index.html").write_text("<!doctype html><title>app</title>")
    (static / "assets" / "app-abc123.js").write_text("console.log(1)")
    (tmp_path / "secret.txt").write_text("nope")
    app = create_app(Settings(db_path=path, static_dir=static, worker_url=None, worker_token=""))
    with TestClient(app) as client:
        assert "<title>app</title>" in client.get("/variants/X12Y").text
        asset = client.get("/assets/app-abc123.js")
        assert asset.status_code == 200 and "immutable" in asset.headers["cache-control"]
        assert "nope" not in client.get("/..%2Fsecret.txt").text
        assert client.get("/api/unknown").status_code == 404
        assert client.post("/api/live/fold", json={"sequence": "MKTAYIAK"}).status_code == 503


def test_client_key_groups_ipv6_by_64():
    from app.live import client_key

    assert client_key("203.0.113.7") == "203.0.113.7"
    assert client_key("2001:db8:1:2:aaaa::1") == client_key("2001:db8:1:2:ffff::9") == "2001:db8:1:2::/64"
    assert client_key("::ffff:203.0.113.7") == "203.0.113.7"
    assert client_key(None) == "unknown"


def test_rate_limiter_forgets_expired_clients_and_refunds():
    from app.live import LiveError, RateLimiter

    rl = RateLimiter(limit=1, window_s=10)
    for i in range(10_001):
        rl.check(f"c{i}", now=0.0)
    rl.check("late", now=100.0)  # every earlier client has expired
    assert len(rl._hits) == 1
    rl.check("x", now=200.0)
    with pytest.raises(LiveError):
        rl.check("x", now=201.0)
    rl.refund("x")
    rl.check("x", now=202.0)


def test_busy_worker_does_not_consume_quota(client_and_worker):
    client, worker = client_and_worker  # live_per_hour=2
    worker.busy = True
    for _ in range(3):
        assert client.post("/api/live/fold", json={"sequence": "MKTAYIAK"}).status_code == 503
    worker.busy = False
    assert client.post("/api/live/fold", json={"sequence": "MKTAYIAK"}).status_code == 200


def test_unicode_digit_search_is_not_a_server_error(client_and_worker):
    client, _ = client_and_worker
    assert client.get("/api/variants", params={"q": "²"}).status_code == 200
