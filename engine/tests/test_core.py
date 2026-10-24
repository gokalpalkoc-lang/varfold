import gzip
import json
import sqlite3

import numpy as np
import pytest
from scipy import stats as sps

from varfold.curate import aggregate_variants, classify, scan_mutations
from varfold.sequences import Mutation, SequenceError, apply_mutations, clean_sequence
from varfold.stats import auc, benjamini_hochberg, empirical_lower_p, roc_curve
from varfold.structure import kabsch, lddt, read_ca, rmsd, superpose_pdb

from synthetic import make_pdb, random_chain


# ----------------------------------------------------------------------------- sequences

def test_mutation_parsing():
    assert Mutation.parse("R856H") == Mutation(856, "R", "H")
    assert Mutation.parse("p.Arg856His") == Mutation(856, "R", "H")
    assert Mutation.parse(" t878a ").id == "T878A"
    for bad in ("R856", "X12A", "R12R", "Arg12Xyz", "RR12H"):
        with pytest.raises(SequenceError):
            Mutation.parse(bad)


def test_multi_letter_codes_are_not_residues():
    with pytest.raises(SequenceError):
        Mutation(5, "DE", "A")


def test_apply_mutations_checks_reference():
    assert apply_mutations("ACDE", [Mutation(11, "C", "W")], offset=10) == "AWDE"
    with pytest.raises(SequenceError, match="reference residue"):
        apply_mutations("ACDE", [Mutation(11, "D", "W")], offset=10)
    with pytest.raises(SequenceError, match="outside"):
        apply_mutations("ACDE", [Mutation(20, "A", "W")], offset=10)


def test_clean_sequence():
    assert clean_sequence(" acd\nEF ") == "ACDEF"
    for bad in ("", "ACDX", ">header", "ACD*"):
        with pytest.raises(SequenceError):
            clean_sequence(bad)
    with pytest.raises(SequenceError, match="limit"):
        clean_sequence("A" * 11, max_length=10)


# ----------------------------------------------------------------------------- curation

@pytest.mark.parametrize("terms,expected", [
    ({"pathogenic"}, "pathogenic"),
    ({"likely pathogenic", "variant of uncertain significance"}, "pathogenic"),
    ({"benign", "likely benign"}, "benign"),
    ({"pathogenic", "likely benign"}, "conflicting"),
    ({"conflicting interpretations of pathogenicity"}, "conflicting"),
    ({"variant of uncertain significance"}, "vus"),
    ({"risk factor"}, "other"),
    (set(), "unlabeled"),
])
def test_classify(terms, expected):
    assert classify(terms) == expected


def _feat(pos, wt, mt, sig=None, sources=("ClinVar",), review="criteria provided, single submitter", **kw):
    f = {"type": "VARIANT", "consequenceType": "missense", "begin": str(pos), "end": str(pos),
         "wildType": wt, "mutatedType": mt, "xrefs": [], **kw}
    if sig:
        f["clinicalSignificances"] = [{"type": sig, "sources": list(sources), "reviewStatus": review}]
    return f


def test_aggregate_variants_groups_and_validates():
    seq = "MKTAYIAK"
    features = [
        _feat(2, "K", "E", "Pathogenic"),
        _feat(2, "K", "E", "Likely pathogenic", review="criteria provided, multiple submitters, no conflicts"),
        _feat(3, "T", "A", "Benign"),
        _feat(3, "T", "A", "Pathogenic", sources=("NCI-TCGA",)),   # somatic-only: ignored
        _feat(4, "Q", "P", "Pathogenic"),                          # reference mismatch (A at 4)
        _feat(5, "Y", "*", "Pathogenic"),                          # nonsense mislabelled as missense
        _feat(99, "A", "G", "Pathogenic"),                         # beyond sequence
        _feat(6, "I", "V"),                                        # no clinical assertion
        {"type": "VARIANT", "consequenceType": "stop gained", "begin": "7", "wildType": "A", "mutatedType": "*"},
    ]
    variants, dropped = aggregate_variants(features, seq)
    by_id = {v["variant_id"]: v for v in variants}
    assert set(by_id) == {"K2E", "T3A", "I6V"}
    assert by_id["K2E"]["label"] == "pathogenic"
    assert by_id["K2E"]["review_stars"] == 2
    assert by_id["K2E"]["assertions"] == "LP:1;P:1"
    assert by_id["T3A"]["label"] == "benign"
    assert by_id["I6V"]["label"] == "unlabeled"
    assert dropped == {"reference_mismatch": 2, "nonsense_or_non_canonical": 1}


def test_scan_mutations():
    assert [m.id for m in scan_mutations("AC", 10, "alanine")] == ["C11A"]
    assert len(scan_mutations("AC", 10, "saturation")) == 38
    assert scan_mutations("AC", 10, "none") == []


# ----------------------------------------------------------------------------- structure

def test_kabsch_recovers_rigid_motion():
    rng = np.random.default_rng(0)
    a = random_chain(50, rng)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))
    b = a @ q.T + np.array([5.0, -3.0, 12.0])
    r, t = kabsch(a, b)
    assert np.allclose(a @ r.T + t, b, atol=1e-8)
    assert rmsd(a, b) < 1e-8
    assert lddt(b, a) == pytest.approx(1.0)


def test_lddt_detects_local_change_rmsd_does_not_mask():
    rng = np.random.default_rng(1)
    a = random_chain(60, rng)
    b = a.copy()
    b[20:26] += rng.normal(0, 3, (6, 3))
    assert lddt(b, a) < 0.97
    assert rmsd(b, a) > 0.3


def test_pdb_roundtrip_and_superpose():
    rng = np.random.default_rng(2)
    seq = "ACDEFGHIKL"
    ca = random_chain(10, rng)
    pdb = make_pdb(seq, ca + 20.0, np.full(10, 77.0))
    model = read_ca(pdb)
    assert list(model.resnum) == list(range(1, 11))
    assert np.allclose(model.bfactor, 77.0)
    moved = read_ca(superpose_pdb(pdb, ca)).coords
    assert np.allclose(moved, ca, atol=2e-3)


# ----------------------------------------------------------------------------- statistics

def test_auc_matches_pairwise_definition():
    rng = np.random.default_rng(3)
    s = np.round(rng.normal(size=60), 1)  # rounding creates ties
    y = rng.random(60) < 0.4
    pairwise = np.mean([(a > b) + 0.5 * (a == b) for a in s[y] for b in s[~y]])
    assert auc(s, y) == pytest.approx(pairwise)
    assert auc(np.array([3, 2, 1, 0]), np.array([1, 1, 0, 0], bool)) == 1.0
    fpr, tpr = roc_curve(s, y)
    assert (fpr[0], tpr[0], fpr[-1], tpr[-1]) == (0, 0, 1, 1)
    assert np.trapezoid(tpr, fpr) == pytest.approx(auc(s, y))


def test_benjamini_hochberg_matches_scipy():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan, 0.001])
    q = benjamini_hochberg(p)
    expected = sps.false_discovery_control(p[~np.isnan(p)], method="bh")
    assert np.allclose(q[~np.isnan(p)], expected)
    assert np.isnan(q[4])


def test_empirical_p():
    null = np.arange(99, dtype=float)
    assert empirical_lower_p(np.array([-1.0]), null)[0] == pytest.approx(1 / 100)
    assert empirical_lower_p(np.array([1000.0]), null)[0] == pytest.approx(1.0)


# ----------------------------------------------------------------------------- pipeline

def test_analysis_and_export_recover_planted_signal(synthetic_run):
    from varfold.analyze import run_analysis
    from varfold.export import export_web_db

    cfg = synthetic_run
    result = run_analysis(cfg)
    disc = {d["metric"]: d for d in result["discrimination"]}
    assert disc["delta_local_plddt"]["auc"] > 0.95
    assert disc["lddt_local"]["auc"] > 0.9
    # Planted distortions are local, so rigid-body motion must not leak into metrics.
    bp = result["breaking_points"]["by_class"]
    assert bp["pathogenic"]["n_breaking_points"] >= 14
    assert bp["benign"]["n_breaking_points"] == 0
    assert (cfg.paths.analysis / "report.md").read_text().startswith("# Toy domain")
    assert (cfg.paths.figures / "fig3_roc.pdf").exists()

    db = export_web_db(cfg)
    con = sqlite3.connect(db)
    n, = con.execute("SELECT COUNT(*) FROM variants").fetchone()
    assert n == 33
    pdb_gz, plddt = con.execute("SELECT pdb_gz, plddt FROM variants WHERE label='benign' LIMIT 1").fetchone()
    wt_gz, = con.execute("SELECT pdb_gz FROM variants WHERE variant_id='WT'").fetchone()
    # Exported variants are superposed onto the WT frame.
    assert rmsd(read_ca(gzip.decompress(pdb_gz).decode()).coords,
                read_ca(gzip.decompress(wt_gz).decode()).coords, superpose=False) < 0.5
    assert len(json.loads(plddt)) == 80
    meta = dict(con.execute("SELECT key, value FROM meta"))
    assert json.loads(meta["construct"])["start"] == 101
    assert con.execute("SELECT COUNT(*) FROM positions").fetchone()[0] == 80


# ----------------------------------------------------------------------------- review regressions

def test_set_bfactors_keeps_full_precision():
    from varfold.folding import _set_bfactors

    pdb = make_pdb("ACD", np.zeros((3, 3)), np.zeros(3))
    out = read_ca(_set_bfactors(pdb, np.array([87.234, 50.5, 99.999])))
    assert np.allclose(out.bfactor, [87.23, 50.5, 100.0])
    with pytest.raises(RuntimeError):
        _set_bfactors(pdb, np.array([1.0, 2.0]))


def test_read_predictions_skips_torn_line(tmp_path):
    from varfold.folding import read_predictions

    p = tmp_path / "predictions.jsonl"
    p.write_text('{"variant_id": "WT", "fingerprint": "a"}\n{"variant_id": "A1C", "fing')
    assert set(read_predictions(p)) == {"WT"}


def test_structure_file_per_fingerprint_and_legacy(tmp_path):
    from varfold.config import RunPaths
    from varfold.folding import fingerprint_dir, structure_file

    paths = RunPaths(tmp_path)
    assert fingerprint_dir("m|recycles=4|fp32") != fingerprint_dir("m|recycles=4|fp16-lm/tf32")
    rec = {"variant_id": "A1C", "structure": "abc/A1C.pdb"}
    assert structure_file(paths, rec) == tmp_path / "structures" / "abc" / "A1C.pdb"
    assert structure_file(paths, {"variant_id": "A1C"}) == tmp_path / "structures" / "A1C.pdb"


def test_changed_construct_is_refused(synthetic_run):
    from varfold.analyze import run_analysis

    synthetic_run.construct.start += 1
    synthetic_run.construct.end += 1
    with pytest.raises(RuntimeError, match="re-run `varfold curate`"):
        run_analysis(synthetic_run)


def test_paired_auc_difference_sanity():
    from varfold.stats import paired_auc_difference

    rng = np.random.default_rng(5)
    y = np.r_[np.ones(80, bool), np.zeros(30, bool)]
    good = y + rng.normal(0, 0.5, len(y))
    same = paired_auc_difference(good, good.copy(), y, n_perm=499)
    assert same["difference"] == 0 and same["p"] == 1.0
    better = paired_auc_difference(good, rng.normal(size=len(y)), y, n_perm=499)
    assert better["difference"] > 0.3 and better["p"] < 0.01 and better["p_delong"] < 0.01
    lo, hi = better["ci95"]
    assert lo < better["difference"] < hi


def test_fisher_odds_ratio_kinds():
    from varfold.stats import fisher_enrichment

    s = np.array([1, 1, 0, 0], bool)
    assert fisher_enrichment(s, np.array([1, 0, 0, 0], bool))["odds_ratio_kind"] == "infinite"
    assert fisher_enrichment(s, np.zeros(4, bool))["odds_ratio_kind"] == "undefined"
    assert fisher_enrichment(np.array([1, 1, 1, 0, 0, 0], bool),
                             np.array([1, 1, 0, 1, 0, 0], bool))["odds_ratio_kind"] == "finite"


def test_compare_groups_reports_raw_medians():
    from varfold.stats import compare_groups

    y = np.array([1, 1, 1, 0, 0, 0], bool)
    res = compare_groups(np.array([-10.0, -8, -9, 0, 1, -1]), y, direction=-1, n_boot=50, seed=0)
    assert res["auc"] == 1.0 and res["median_pos"] == -9.0 and res["median_neg"] == 0.0
