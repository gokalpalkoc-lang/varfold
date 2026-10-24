import type { Meta } from "../lib/api";
import { date, num } from "../lib/format";

const STEPS = [
  { n: "1", title: "Curate", body: "Canonical sequence from UniProtKB; missense variants with clinical assertions from the EBI Proteins API (ClinVar, UniProt curation, Ensembl), checked residue-by-residue against the reference." },
  { n: "2", title: "Fold", body: "Every sequence folded with ESMFold on a local GPU; atomic coordinates and per-residue pLDDT/pTM logged with the exact model settings." },
  { n: "3", title: "Validate", body: "The wild-type prediction is benchmarked against an experimental crystal structure, with residues matched through SIFTS." },
  { n: "4", title: "Compare", body: "Each variant is compared with the wild type: local and global ΔpLDDT, Cα RMSD and superposition-free lDDT-Cα." },
  { n: "5", title: "Test", body: "Breaking points are called against an empirical noise floor; pathogenic vs benign discrimination is tested against ESM-2 and BLOSUM62 baselines." },
];

export function Methods({ meta }: { meta: Meta }) {
  const s = meta.stats;
  const nf = s.noise_floor;
  const bp = s.breaking_points;
  const v = meta.validation;
  const filters = meta.provenance.filters as { include_classes?: string[]; min_review_stars?: number; scan?: string };
  return (
    <div className="stack">
      <section className="stack-sm">
        <h1>Methods</h1>
        <p className="lede">How the dataset behind this site was built, what each metric means, and where its limits are.</p>
      </section>

      <section className="grid-tiles" aria-label="Pipeline">
        {STEPS.map((st) => (
          <div className="tile" key={st.n}>
            <div className="eyebrow">Step {st.n}</div>
            <h3 style={{ margin: "4px 0 6px" }}>{st.title}</h3>
            <p className="small" style={{ color: "var(--ink-2)" }}>{st.body}</p>
          </div>
        ))}
      </section>

      <article className="card prose" style={{ maxWidth: "none" }}>
        <h2 style={{ marginTop: 0 }}>Target and dataset</h2>
        <p>
          {meta.target.title} ({meta.target.gene}, UniProt <a href={`https://www.uniprot.org/uniprotkb/${meta.target.uniprot}`} target="_blank" rel="noreferrer">{meta.target.uniprot}</a>
          {meta.provenance.uniprot_entry_version ? `, entry version ${meta.provenance.uniprot_entry_version}` : ""}), residues{" "}
          {meta.construct.start}–{meta.construct.end} ({meta.construct.length} aa). All numbering is UniProt canonical; older
          literature on this protein often uses a shifted numbering.
        </p>
        <p>
          Sources were retrieved on {date(meta.provenance.retrieved_at)} and frozen as raw snapshots; their SHA-256 digests are
          listed below. Records whose stated reference residue does not match the canonical sequence (typically mapped through
          another isoform) were dropped, not re-labelled. Assertions for one substitution are merged as follows: any conflicting
          assertion, or both pathogenic and benign, gives <i>conflicting</i>; otherwise P/LP gives <i>pathogenic</i> and B/LB
          gives <i>benign</i> (VUS assertions alongside are tolerated); VUS-only gives <i>VUS</i>. Somatic-only assertions are ignored.
        </p>
        <ul>
          <li>Classes included: {(filters.include_classes ?? []).join(", ") || "—"}; minimum ClinVar review: {filters.min_review_stars ?? 0}★; scan: {filters.scan ?? "none"}.</li>
          <li>Dropped during curation: {Object.entries(meta.provenance.dropped).map(([k, n]) => `${k.replace(/_/g, " ")} ${n}`).join(" · ") || "none"}.</li>
          {Object.entries(meta.provenance.sources).map(([name, src]) => (
            <li key={name}><span className="mono">{src.file}</span> — <a href={src.url} target="_blank" rel="noreferrer">{name}</a>, sha256 <span className="mono">{src.sha256.slice(0, 16)}…</span></li>
          ))}
        </ul>

        <h2>Structure prediction</h2>
        <p>
          Model <span className="mono">{meta.model.folding}</span> with {meta.model.num_recycles} recycles
          {meta.model.precision ? ` (${meta.model.precision === "fp16-lm" ? "ESM-2 language model in fp16, folding trunk in fp32" : "fp32"})` : ""}
          {meta.model.device ? ` on ${meta.model.device}` : ""}. Settings fingerprint: <span className="mono">{meta.model.fingerprint}</span>.
          Every variant is folded as the full construct with a single substitution and compared only with a wild type folded under
          identical settings. Downloadable structures are rigidly superposed onto the wild-type frame (Cα Kabsch).
        </p>
        {v && (
          <p>
            Validation: against {v.pdb_id} chain {v.chain} ({v.n_aligned} residues, {Math.round(v.coverage * 100)}% of the construct)
            the wild-type prediction has Cα RMSD {num(v.ca_rmsd, 2)} Å and lDDT-Cα {num(v.lddt_ca, 3)}; per-residue pLDDT correlates
            with per-residue lDDT at Spearman ρ = {num(v.plddt_vs_lddt_spearman.rho, 2)}.
          </p>
        )}

        <h2>Metrics</h2>
        <ul>
          <li><b>Δ local pLDDT</b> — mean change in per-residue (Cα) pLDDT within ±{bp.local_window} residues of the site. The primary metric.</li>
          <li><b>Δ mean pLDDT / Δ site pLDDT / Δ pTM</b> — whole-domain, single-residue and global-fold confidence changes.</li>
          <li><b>Cα RMSD</b> — after optimal superposition on the wild-type prediction.</li>
          <li><b>lDDT-Cα</b> — superposition-free agreement of all Cα–Cα distances under 15 Å (thresholds 0.5, 1, 2, 4 Å); reported globally and within the local window.</li>
          <li><b>ESM-2 LLR</b> — log p(mutant) − log p(wild type) at the masked site under <span className="mono">{meta.model.llr}</span>: the zero-shot, sequence-only baseline.</li>
          <li><b>BLOSUM62</b> — substitution-matrix score of the amino-acid swap: the trivial chemistry baseline. Benign variants are dominated by conservative swaps, so a structural metric is only interesting if it beats this.</li>
          <li><b>Site priors</b> — wild-type pLDDT and Cα burial at the site. They ignore which residue is substituted, so any structural metric should beat them.</li>
        </ul>

        <h2>Breaking points</h2>
        <p>
          ESMFold is deterministic, so a noise floor cannot come from replicate runs. Instead, for every variant the same windowed
          ΔpLDDT statistic is evaluated at windows lying entirely more than 15 Å from the mutation (n = {nf.n.toLocaleString()} windows), where a
          substitution should have no direct effect. A variant's local ΔpLDDT gets an empirical one-sided p-value against that
          distribution, corrected with Benjamini–Hochberg; variants with q &lt; {bp.fdr} are <b>breaking points</b>. The noise floor's 1st
          percentile is {num(nf.p01, 2)} pLDDT points.
        </p>

        <h2>Statistics</h2>
        <p>
          Pathogenic vs benign discrimination is summarised by ROC AUC (Mann–Whitney identity) with 95% CIs from a stratified
          bootstrap, a two-sided Mann–Whitney U test, Cliff's δ, and Benjamini–Hochberg q-values across metrics. Each structural
          metric is compared with the ESM-2 and BLOSUM62 baselines and the site prior on the same variants: the AUC difference
          gets a DeLong 95% CI, and its p-value comes from a paired permutation test, because DeLong's test is liberal with a
          benign class this small. Whether pathogenic variants become breaking points more often than benign ones is a one-sided
          Fisher exact test at the variant level. Analysis generated {date(s.generated_at)}.
        </p>

        <h2>Limitations</h2>
        <ul>
          <li>ESMFold was trained to predict structures, not mutation effects; single substitutions often leave predictions unchanged even when they are pathogenic.</li>
          <li>pLDDT measures model confidence, not thermodynamic stability or function. Pathogenic mechanisms such as loss of ligand or co-activator binding need not perturb the fold.</li>
          <li>The benign class is small, so AUC estimates are uncertain; clinical databases also carry ascertainment bias toward well-studied regions.</li>
          <li>Research use only. Nothing here is a clinical interpretation.</li>
        </ul>

        <h2>Programmatic access</h2>
        <p>
          Everything on this site is available through a read-only JSON API documented at <a href="/api/docs">/api/docs</a>, for example{" "}
          <span className="mono">GET /api/variants?label=pathogenic&amp;breaking=true</span> or{" "}
          <span className="mono">GET /api/variants/{"{id}"}/structure.pdb</span>.
        </p>

        <h2>References</h2>
        <ul className="small">
          <li>Lin Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. <i>Science</i> 379, 1123–1130 (2023).</li>
          <li>Meier J. et al. Language models enable zero-shot prediction of the effects of mutations on protein function. <i>NeurIPS</i> (2021).</li>
          <li>Henikoff S. &amp; Henikoff J.G. Amino acid substitution matrices from protein blocks. <i>PNAS</i> 89, 10915–10919 (1992).</li>
          <li>Mariani V. et al. lDDT: a local superposition-free score for comparing protein structures and models. <i>Bioinformatics</i> 29, 2722–2728 (2013).</li>
          <li>Landrum M.J. et al. ClinVar: improvements to accessing data. <i>Nucleic Acids Res.</i> 48, D835–D844 (2020).</li>
          <li>The UniProt Consortium. UniProt: the Universal Protein Knowledgebase. <i>Nucleic Acids Res.</i> (2025).</li>
          <li>Rego N. &amp; Koes D. 3Dmol.js: molecular visualization with WebGL. <i>Bioinformatics</i> 31, 1322–1324 (2015).</li>
        </ul>
      </article>
    </div>
  );
}
