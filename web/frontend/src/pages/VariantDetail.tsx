import { useMemo, useState } from "react";
import { Link, useParams } from "react-router";
import { api, useAsync, type Meta, type VariantDetail as Detail } from "../lib/api";
import { labelVar } from "../lib/colors";
import { num, pval } from "../lib/format";
import { MolViewer } from "../components/MolViewer";
import { ResidueChart } from "../components/ResidueChart";
import { BreakingBadge, ClassBadge, DeltaLegend, ErrorNote, LineLegend, Loading, PlddtLegend, StatTile } from "../components/ui";

const THREE: Record<string, string> = {
  A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile",
  L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val",
};

export function VariantDetail({ meta }: { meta: Meta }) {
  const { id = "" } = useParams();
  const variant = useAsync(() => api.variant(id), id);
  const pdb = useAsync(() => api.structure(id), `pdb-${id}`);
  const wtPdb = useAsync(() => api.structure("WT"), "WT");

  if (variant.error) {
    return (
      <div className="stack">
        <ErrorNote error={variant.error} />
        <Link to="/variants">← Back to all variants</Link>
      </div>
    );
  }
  if (!variant.data || variant.data.variant_id !== id) return <Loading what={id} />;
  return (
    <>
      {pdb.error && <ErrorNote error={pdb.error} />}
      <DetailBody v={variant.data} meta={meta} pdb={pdb.data} wtPdb={wtPdb.data} />
    </>
  );
}

function DetailBody({ v, meta, pdb, wtPdb }: { v: Detail; meta: Meta; pdb?: string; wtPdb?: string }) {
  const [colorBy, setColorBy] = useState<"plddt" | "delta">("plddt");
  const [ghost, setGhost] = useState(true);
  const [hover, setHover] = useState<number | null>(null);
  const start = meta.construct.start;
  const w = meta.stats.breaking_points.local_window;
  const wt = meta.wild_type.plddt;
  const delta = useMemo(() => v.plddt.map((p, i) => p - wt[i]), [v.plddt, wt]);
  const isWt = v.variant_id === "WT";
  const resi = v.position !== null ? v.position - start + 1 : undefined;
  const span = Math.max(5, ...delta.map(Math.abs));

  return (
    <div className="stack">
      <section className="stack-sm">
        <Link to="/variants" className="small">← All variants</Link>
        <div className="row" style={{ gap: 14 }}>
          <h1 style={{ fontFamily: "var(--mono)" }}>{v.variant_id}</h1>
          {!isWt && <span className="muted">p.{THREE[v.wt!]}{v.position}{THREE[v.mt!]}</span>}
          <ClassBadge label={v.label} />
          {v.breaking_point ? <BreakingBadge q={v.local_q} /> : null}
        </div>
        {v.diseases && <p className="muted">{v.diseases.split(";").join(" · ")}</p>}
        <ClinicalLinks v={v} />
      </section>

      <section className="grid-2">
        <div className="card">
          <div className="card-head">
            <h2>Predicted structure</h2>
            <div className="row">
              <div className="segmented" role="group" aria-label="Colour by">
                <button aria-pressed={colorBy === "plddt"} onClick={() => setColorBy("plddt")}>pLDDT</button>
                <button aria-pressed={colorBy === "delta"} onClick={() => setColorBy("delta")} disabled={isWt}>Δ vs WT</button>
              </div>
              {!isWt && (
                <label className="toggle"><input type="checkbox" checked={ghost} onChange={(e) => setGhost(e.target.checked)} />WT ghost</label>
              )}
            </div>
          </div>
          <MolViewer
            pdb={pdb} values={colorBy === "delta" ? delta : v.plddt} colorBy={colorBy}
            overlayPdb={!isWt && ghost ? wtPdb : undefined} highlightResi={resi}
            highlightLabel={isWt ? undefined : v.variant_id} numberOffset={start - 1} sequence={v.sequence}
            label={`Predicted structure of ${v.variant_id}${!isWt && ghost ? " superposed on the wild type" : ""}`}
          />
          <div style={{ marginTop: 12 }}>{colorBy === "plddt" ? <PlddtLegend /> : <DeltaLegend />}</div>
          <div className="row small" style={{ marginTop: 10 }}>
            <a href={`/api/variants/${v.variant_id}/structure.pdb`} download={`${v.variant_id}.pdb`}>Download PDB</a>
            <span className="muted">· superposed on the wild-type frame · B-factor = pLDDT</span>
          </div>
        </div>

        <div className="stack-sm">
          <div className="grid-tiles" style={{ gridTemplateColumns: "repeat(2, minmax(0, 1fr))" }}>
            <StatTile label={`Δ local pLDDT (±${w})`} value={num(v.delta_local_plddt, 2, true)}
                      sub={v.local_q != null ? `noise-floor q = ${pval(v.local_q)}` : undefined} />
            <StatTile label="Δ mean pLDDT" value={num(v.delta_mean_plddt, 2, true)} sub={`mutant ${num(v.mean_plddt)}`} />
            <StatTile label="Site pLDDT" value={num(v.site_plddt)} sub={`WT ${num(v.wt_site_plddt)}`} />
            <StatTile label="Cα RMSD to WT" value={v.rmsd_ca != null ? `${num(v.rmsd_ca, 2)} Å` : "—"} sub={`local lDDT ${num(v.lddt_local, 3)}`} />
            <StatTile label="ESM-2 LLR" value={num(v.llr, 2, true)}
                      sub={`sequence baseline · BLOSUM62 ${v.blosum62 ?? "—"}`} />
            <StatTile label="pTM" value={num(v.ptm, 3)} sub={v.delta_ptm != null ? `Δ ${num(v.delta_ptm, 3, true)}` : undefined} />
          </div>
          {!isWt && (
            <p className="callout">
              {v.breaking_point
                ? <>ESMFold loses local confidence around residue {v.position} beyond the noise floor measured &gt;15 Å from each
                    mutation (FDR {meta.stats.breaking_points.fdr}). This marks a <b>structural breaking point</b> for the model.</>
                : <>Local confidence around residue {v.position} stays within the noise floor: ESMFold predicts essentially the
                    same fold. pLDDT is model confidence, not thermodynamic stability.</>}
            </p>
          )}
          {v.siblings && v.siblings.length > 0 && (
            <div className="card" style={{ padding: 16 }}>
              <h3 style={{ marginBottom: 10 }}>Other substitutions at {v.wt}{v.position}</h3>
              <div className="row">
                {v.siblings.map((s) => (
                  <Link key={s.variant_id} to={`/variants/${s.variant_id}`} className="chip" style={{ textDecoration: "none" }}
                        title={`Δ local pLDDT ${num(s.delta_local_plddt, 2, true)}`}>
                    <span className="dot" style={{ background: labelVar(s.label) }} />
                    <span className="mono">{s.variant_id}</span>
                  </Link>
                ))}
              </div>
            </div>
          )}
        </div>
      </section>

      {!isWt && (
        <section className="card">
          <div className="card-head">
            <h2>Per-residue confidence</h2>
            <LineLegend items={[{ label: "Wild type", color: "var(--wt-line)" }, { label: v.variant_id, color: "var(--mutant)" }]} />
          </div>
          <ResidueChart
            title={`pLDDT per residue, wild type vs ${v.variant_id}`} start={start} sequence={meta.construct.sequence}
            lines={[
              { key: "wt", label: "Wild type", values: wt, color: "var(--wt-line)" },
              { key: "mut", label: v.variant_id, values: v.plddt, color: "var(--mutant)" },
            ]}
            yDomain={[0, 100]} yTicks={[0, 50, 70, 90, 100]} yLabel="pLDDT" height={210}
            marker={v.position ?? undefined} windowRange={v.position ? [v.position - w, v.position + w] : undefined}
            hover={hover} onHover={setHover}
          />
          <ResidueChart
            title={`pLDDT change per residue for ${v.variant_id}`} start={start} sequence={meta.construct.sequence}
            stems={{ label: "ΔpLDDT", values: delta, color: (d) => (d < 0 ? "var(--neg)" : "var(--pos)") }}
            yDomain={[-Math.ceil(span), Math.ceil(span)]} yLabel="ΔpLDDT" height={170}
            marker={v.position ?? undefined} windowRange={v.position ? [v.position - w, v.position + w] : undefined}
            hover={hover} onHover={setHover} format={(x) => (x > 0 ? `+${x.toFixed(1)}` : x.toFixed(1))}
          />
          <p className="small muted" style={{ marginTop: 8 }}>
            Shaded band: the ±{w}-residue window used for Δ local pLDDT. Largest single-residue drop: {num(v.min_delta_plddt, 1)} at
            residue {v.min_delta_position ?? "—"}; {v.n_residues_drop5 ?? 0} residues lose more than 5 points.
          </p>
          <details style={{ marginTop: 10 }}>
            <summary>Show as table</summary>
            <div className="table-wrap" style={{ maxHeight: 360 }}>
              <table>
                <thead><tr><th>Residue</th><th className="num">WT pLDDT</th><th className="num">{v.variant_id} pLDDT</th><th className="num">Δ</th></tr></thead>
                <tbody>
                  {v.plddt.map((p, i) => (
                    <tr key={i}>
                      <td className="mono">{v.sequence[i]}{start + i}</td><td className="num">{num(wt[i])}</td>
                      <td className="num">{num(p)}</td><td className="num">{num(delta[i], 1, true)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </details>
        </section>
      )}
    </div>
  );
}

function ClinicalLinks({ v }: { v: Detail }) {
  const clinvar = (v.clinvar_ids ?? "").split(";").filter(Boolean);
  const dbsnp = (v.dbsnp_ids ?? "").split(";").filter(Boolean);
  if (!clinvar.length && !dbsnp.length && !v.assertions) return null;
  const stars = v.review_stars != null ? "★".repeat(v.review_stars) + "☆".repeat(Math.max(0, 4 - v.review_stars)) : null;
  return (
    <div className="row small" style={{ color: "var(--ink-2)" }}>
      {v.assertions && <span>Assertions: <span className="mono">{v.assertions.replace(/;/g, " · ")}</span></span>}
      {stars && <span title="Best ClinVar review status">ClinVar review {stars}</span>}
      {v.uniprot_curated ? <span>UniProt-curated</span> : null}
      {clinvar.slice(0, 4).map((id) => (
        <a key={id} href={`https://www.ncbi.nlm.nih.gov/clinvar/${id}/`} target="_blank" rel="noreferrer">{id}</a>
      ))}
      {clinvar.length > 4 && <span className="muted">+{clinvar.length - 4} more</span>}
      {dbsnp.slice(0, 2).map((id) => (
        <a key={id} href={`https://www.ncbi.nlm.nih.gov/snp/${id}`} target="_blank" rel="noreferrer">{id}</a>
      ))}
    </div>
  );
}
