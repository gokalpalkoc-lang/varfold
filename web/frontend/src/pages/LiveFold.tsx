import { useEffect, useState, type FormEvent } from "react";
import { Link } from "react-router";
import { api, ApiError, type LiveResult, type LiveStatus, type Meta } from "../lib/api";
import { num } from "../lib/format";
import { MolViewer } from "../components/MolViewer";
import { ResidueChart } from "../components/ResidueChart";
import { LineLegend, PlddtLegend, StatTile } from "../components/ui";

const AA = /^[ACDEFGHIKLMNPQRSTVWY]+$/;

function useLiveStatus(): LiveStatus | undefined {
  const [status, setStatus] = useState<LiveStatus>();
  useEffect(() => {
    let alive = true;
    const poll = () => api.liveStatus().then((s) => alive && setStatus(s), () => alive && setStatus({ available: false, configured: true }));
    poll();
    const t = window.setInterval(poll, 20_000);
    return () => {
      alive = false;
      window.clearInterval(t);
    };
  }, []);
  return status;
}

export function LiveFold({ meta }: { meta: Meta }) {
  const status = useLiveStatus();
  const { start, sequence: wtSeq } = meta.construct;
  const [seq, setSeq] = useState(wtSeq);
  const [mutation, setMutation] = useState("");
  const [mutationError, setMutationError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [elapsed, setElapsed] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<LiveResult | null>(null);
  const [wtPdb, setWtPdb] = useState<string>();

  useEffect(() => {
    if (!busy) return;
    const t0 = Date.now();
    const t = window.setInterval(() => setElapsed((Date.now() - t0) / 1000), 250);
    return () => window.clearInterval(t);
  }, [busy]);
  useEffect(() => {
    api.structure("WT").then(setWtPdb, () => undefined);
  }, []);

  const clean = seq.replace(/\s+/g, "").toUpperCase();
  const maxLen = status?.max_length ?? meta.live.max_length;
  const invalid = clean.length === 0 ? "Enter a sequence." : !AA.test(clean) ? "Only the 20 standard amino-acid letters are accepted."
    : clean.length > maxLen ? `The live engine accepts up to ${maxLen} residues.` : null;

  const applyMutation = () => {
    const m = /^([A-Z])(\d+)([A-Z])$/.exec(mutation.trim().toUpperCase());
    if (!m) return setMutationError("Use the form R856H (UniProt numbering).");
    const [, from, posStr, to] = m;
    const i = Number(posStr) - start;
    if (clean.length !== wtSeq.length) return setMutationError("Quick mutate works on the reference-length construct.");
    if (i < 0 || i >= clean.length) return setMutationError(`Residue must be within ${start}–${start + clean.length - 1}.`);
    if (clean[i] !== from) return setMutationError(`Residue ${posStr} is ${clean[i]}, not ${from}.`);
    if (!AA.test(to)) return setMutationError("Unknown target residue.");
    setSeq(clean.slice(0, i) + to + clean.slice(i + 1));
    setMutation("");
    setMutationError(null);
  };

  const submit = async (e: FormEvent) => {
    e.preventDefault();
    if (invalid) return;
    setBusy(true);
    setError(null);
    setResult(null);
    setElapsed(0);
    try {
      setResult(await api.liveFold(clean));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Request failed — check your connection.");
    } finally {
      setBusy(false);
    }
  };

  const diff = clean.length === wtSeq.length ? [...clean].flatMap((c, i) => (c !== wtSeq[i] ? [`${wtSeq[i]}${start + i}${c}`] : [])) : null;
  const online = status?.available;

  return (
    <div className="stack">
      <section className="stack-sm">
        <h1>Live fold</h1>
        <p className="lede">
          Fold any sequence on demand. Requests are routed through a private Tailscale tunnel to the research
          workstation's GPU; sequences already in the dataset are answered instantly from precomputed results.
        </p>
        <div className="row">
          <span className={`status ${status === undefined ? "" : online ? "online" : "offline"}`} role="status">
            {status === undefined ? "Checking live engine…"
              : online ? `Live engine online${status.busy ? " · busy" : ""}${status.device ? ` · ${status.device}` : ""}`
              : status.configured === false ? "Live folding is not enabled on this deployment"
              : "Live engine offline — precomputed variants still work"}
          </span>
        </div>
      </section>

      <form className="card stack-sm" onSubmit={submit}>
        <div className="card-head" style={{ marginBottom: 0 }}>
          <h2>Sequence</h2>
          <span className="small muted">{clean.length} / {maxLen} residues{diff && diff.length ? ` · ${diff.length} change${diff.length > 1 ? "s" : ""} vs WT` : diff ? " · identical to WT" : ""}</span>
        </div>
        <label className="sr-only" htmlFor="seq">Protein sequence</label>
        <textarea id="seq" rows={6} value={seq} onChange={(e) => setSeq(e.target.value)} spellCheck={false} aria-invalid={!!invalid} />
        <div className="row">
          <label className="sr-only" htmlFor="mut">Point mutation</label>
          <input id="mut" type="text" placeholder="Quick mutate, e.g. R856H" value={mutation} maxLength={10}
                 onChange={(e) => setMutation(e.target.value)} onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); applyMutation(); } }} />
          <button type="button" className="btn" onClick={applyMutation}>Apply</button>
          <button type="button" className="btn" onClick={() => { setSeq(wtSeq); setMutationError(null); }}>Reset to WT</button>
          <span style={{ flex: 1 }} />
          <button type="submit" className="btn primary" disabled={!!invalid || busy}>
            {busy ? `Folding… ${elapsed.toFixed(0)} s` : "Fold structure"}
          </button>
        </div>
        {(mutationError || (invalid && clean.length > 0)) && <p className="small error" role="alert">{mutationError ?? invalid}</p>}
        {diff && diff.length > 0 && <p className="small muted">Changes vs WT: <span className="mono">{diff.slice(0, 12).join(", ")}{diff.length > 12 ? " …" : ""}</span></p>}
      </form>

      {error && <p className="callout warn" role="alert">{error}</p>}
      {result?.source === "dataset" && (
        <p className="callout">
          This sequence is part of the precomputed dataset ({result.variant_id === "WT" ? "the wild type" : result.variant_id}).{" "}
          <Link to={result.variant_id === "WT" ? "/" : `/variants/${result.variant_id}`}>Open the full analysis →</Link>
        </p>
      )}
      {result?.source === "live" && <LiveResultView r={result} meta={meta} wtPdb={wtPdb} />}
    </div>
  );
}

function LiveResultView({ r, meta, wtPdb }: { r: Extract<LiveResult, { source: "live" }>; meta: Meta; wtPdb?: string }) {
  const comparable = r.mutations !== null;
  const start = meta.construct.start;
  const wt = meta.wild_type.plddt;
  const site = r.mutations?.length === 1 ? Number(r.mutations[0].slice(1, -1)) : undefined;
  return (
    <section className="stack">
      <div className="grid-tiles">
        <StatTile label="Mean pLDDT" value={num(r.mean_plddt)} sub={comparable ? `WT ${num(meta.stats.wild_type.mean_plddt)}` : undefined} />
        <StatTile label="pTM" value={num(r.ptm, 3)} />
        <StatTile label="Inference time" value={`${num(r.runtime_s, 1)} s`} sub={r.cached ? "served from cache" : r.device ?? undefined} />
        <StatTile label="Length" value={r.sequence.length} sub={r.mutations?.length ? r.mutations.join(", ").slice(0, 40) : undefined} />
      </div>
      {comparable && !r.comparable && (
        <p className="callout warn">
          The live engine is running with different settings from the dataset ({meta.model.fingerprint}), so small
          differences from the wild-type values may come from the settings rather than the sequence.
        </p>
      )}
      <div className="card">
        <div className="card-head"><h2>Predicted structure</h2>{comparable && <span className="small muted">grey ghost: wild type</span>}</div>
        <MolViewer pdb={r.pdb} overlayPdb={comparable ? wtPdb : undefined}
                   highlightResi={site !== undefined ? site - start + 1 : undefined}
                   highlightLabel={site !== undefined ? r.mutations![0] : undefined}
                   numberOffset={comparable ? start - 1 : 0} sequence={r.sequence} label="Live ESMFold prediction coloured by pLDDT" />
        <div style={{ marginTop: 12 }}><PlddtLegend /></div>
      </div>
      <div className="card">
        <div className="card-head">
          <h2>Per-residue confidence</h2>
          {comparable && <LineLegend items={[{ label: "Wild type", color: "var(--wt-line)" }, { label: "This sequence", color: "var(--mutant)" }]} />}
        </div>
        <ResidueChart
          title="pLDDT per residue for the live prediction" start={comparable ? start : 1} sequence={r.sequence}
          lines={[
            ...(comparable ? [{ key: "wt", label: "Wild type", values: wt, color: "var(--wt-line)" }] : []),
            { key: "live", label: comparable ? "This sequence" : "pLDDT", values: r.plddt, color: "var(--mutant)" },
          ]}
          yDomain={[0, 100]} yTicks={[0, 50, 70, 90, 100]} yLabel="pLDDT" marker={site}
        />
      </div>
    </section>
  );
}
