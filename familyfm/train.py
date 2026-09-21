"""Fit and evaluate the target-preference comparator.

Split is compound-disjoint: a ligand's InChIKey decides which side it lands on, so
no ligand appears in both fit and test. Ties are already dropped at pairing. Every
comparison is entered twice with the label inverted, and every prediction averages
both target orders, so compare(A,B) + compare(B,A) is 1 and label balance is
exactly 0.5000.

A pairs file from `familyfm.pairs_family_mode --mode both` carries a `pair_type`
column, `cross` for two targets from different families and `same` for two
different targets sharing one. Both kinds are fitted together and measured
SEPARATELY, in MEASURED_PERFORMANCE_BY_PAIR_TYPE. The pooled block is kept for
provenance only: the holdout is mostly same-family, so a pooled accuracy can rise
while the cross-family one falls, and it must never be quoted.

    python -m familyfm.train --pairs data/pairs/all_both.csv \
        --out models/xfam_v2_samefamily
"""
import argparse
import hashlib
import json
import os
import sys
import time

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import embed as emb                                            # noqa: E402
from features import ligand_features, sls_rows, SEQ_DIMS       # noqa: E402

BANDS = (0.50, 0.60, 0.70, 0.80, 0.90)


def test_side(inchikey, fraction, salt="familyfm"):
    h = hashlib.blake2b((salt + str(inchikey)).encode(), digest_size=8).digest()
    return (int.from_bytes(h, "big") % 10000) < int(fraction * 10000)


def cumulative(strength, correct):
    out = {}
    for c in BANDS:
        m = strength >= c
        out["%.2f" % c] = {
            "coverage": round(float(m.mean()), 4),
            "accuracy": None if not m.any() else round(float(correct[m].mean()), 4),
            "n": int(m.sum())}
    return out


def measure(g, per_pair=True, per_family=False):
    """Held-out figures for one set of scored comparisons."""
    out = {"held_out_comparisons": int(len(g)),
           "held_out_ligands": int(g.inchikey.nunique()),
           "accuracy": round(float(g.correct.mean()), 4),
           "strength_cumulative": cumulative(g.strength.values, g.correct.values),
           "by_endpoint": {k: {"n": int(len(h)),
                               "accuracy": round(float(h.correct.mean()), 4)}
                           for k, h in g.groupby("assay_measure")}}
    fa, fb = g.family_a.str.split("|").str[0], g.family_b.str.split("|").str[0]
    if per_pair:
        # key on the SORTED pair, so both orientations of a pairing count together
        key = np.array([" | ".join(sorted(k)) for k in zip(fa, fb)])
        pp = {}
        for kk, h in g.groupby(key):
            if len(h) >= 30:
                pp[kk] = {"n": int(len(h)), "accuracy": round(float(h.correct.mean()), 4)}
        out["by_family_pair"] = pp
    if per_family:
        # a same-family comparison is filed under every family its two targets
        # share (a Reader|Writer pair counts under both); `targets` is how many
        # distinct proteins that figure rests on
        pf = {}
        shared = [set(x.split("|")) & set(y.split("|"))
                  for x, y in zip(g.family_a, g.family_b)]
        for fam in sorted(set().union(*shared)):
            h = g[[fam in s for s in shared]]
            pf[fam] = {"n": int(len(h)), "accuracy": round(float(h.correct.mean()), 4),
                       "targets": int(len(set(h.uniprot_a) | set(h.uniprot_b))),
                       "scored": bool(len(h) >= 200)}
        out["by_family"] = pf
    gap = np.abs(g.pic50_a - g.pic50_b)
    gs = {}
    for lo, hi, nm in [(0, .5, "under 0.5 log"), (.5, 1, "0.5 to 1 log"),
                       (1, 2, "1 to 2 logs"), (2, 99, "over 2 logs")]:
        m = (gap >= lo) & (gap < hi)
        if m.sum():
            gs[nm] = {"n": int(m.sum()), "accuracy": round(float(g.correct[m].mean()), 4)}
    out["gap_strata"] = gs
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--trees", type=int, default=300)
    ap.add_argument("--min-samples-leaf", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    d = pd.read_csv(a.pairs, low_memory=False)
    d["is_test"] = [test_side(k, a.test_fraction) for k in d.inchikey]
    fit, tst = d[~d.is_test].copy(), d[d.is_test].copy()
    print("comparisons %d: fit %d, test %d" % (len(d), len(fit), len(tst)), flush=True)
    print("ligands: fit %d, test %d, shared %d"
          % (fit.inchikey.nunique(), tst.inchikey.nunique(),
             len(set(fit.inchikey) & set(tst.inchikey))), flush=True)

    # ---- features ----
    seqs = sorted(set(d.sequence_a) | set(d.sequence_b))
    print("embedding %d sequences with %s" % (len(seqs), emb.CHECKPOINT), flush=True)
    t0 = time.time()
    SV = {}
    for i, s in enumerate(seqs, 1):
        SV[s] = np.asarray(emb.embed_sequence(s), dtype=np.float32)[:SEQ_DIMS]
        if i % 250 == 0:
            print("  %d/%d" % (i, len(seqs)), flush=True)
    print("embedded in %.0f s" % (time.time() - t0), flush=True)

    ligs = sorted(set(d.smiles))
    print("featurising %d ligands" % len(ligs), flush=True)
    LF = ligand_features(ligs)
    li = {s: i for i, s in enumerate(ligs)}

    def matrix(frame):
        sa = np.stack([SV[s] for s in frame.sequence_a])
        sb = np.stack([SV[s] for s in frame.sequence_b])
        lg = LF[[li[s] for s in frame.smiles]]
        return sa, lg, sb

    fa, fl, fb = matrix(fit)
    ta, tl, tb = matrix(tst)
    y = (fit.pic50_a.values > fit.pic50_b.values).astype(int)

    # swap augmentation: both orders, label inverted
    X = np.vstack([sls_rows(fa, fl, fb), sls_rows(fb, fl, fa)])
    Y = np.concatenate([y, 1 - y])
    print("fit matrix %s, label balance %.4f" % (X.shape, Y.mean()), flush=True)

    rf = RandomForestClassifier(n_estimators=a.trees,
                                min_samples_leaf=a.min_samples_leaf,
                                max_features="sqrt", random_state=a.seed,
                                n_jobs=-1)
    t0 = time.time()
    rf.fit(X, Y)
    fit_s = time.time() - t0
    print("fitted in %.0f s" % fit_s, flush=True)
    del X, Y

    c1 = int(np.where(rf.classes_ == 1)[0][0])

    def score(rfm, sa, lg, sb):
        p1 = rfm.predict_proba(sls_rows(sa, lg, sb))[:, c1]
        p2 = rfm.predict_proba(sls_rows(sb, lg, sa))[:, c1]
        return 0.5 * (p1 + (1.0 - p2))

    P = score(rf, ta, tl, tb)
    lab = (tst.pic50_a.values > tst.pic50_b.values)
    correct = (P > 0.5) == lab
    strength = np.maximum(P, 1 - P)
    acc = float(correct.mean())
    print("\nMODEL held-out accuracy %.4f on %d comparisons%s"
          % (acc, len(tst), ", pooled, do not quote" if "pair_type" in tst.columns
             else ""), flush=True)

    tst = tst.assign(p=P, correct=correct, strength=strength)
    by_type = None
    if "pair_type" in tst.columns:
        by_type = {}
        for kind in ("cross", "same"):
            g = tst[tst.pair_type == kind]
            if len(g):
                by_type[kind] = measure(g, per_pair=(kind == "cross"),
                                        per_family=(kind == "same"))
                print("  %-5s arm accuracy %.4f on %d comparisons"
                      % (kind, by_type[kind]["accuracy"], len(g)), flush=True)
    pooled = measure(tst, per_pair=True, per_family=False)
    per_pair, gap_strata = pooled["by_family_pair"], pooled["gap_strata"]

    joblib.dump(rf, os.path.join(a.out, "familyfm_selectivity.joblib"), compress=3)
    idx = {emb.sequence_key(s): i for i, s in enumerate(seqs)}
    np.savez_compressed(os.path.join(a.out, "sequence_vectors.npz"),
                        vectors=np.stack([SV[s] for s in seqs]))
    json.dump({"seq_key_to_row": idx}, open(os.path.join(a.out, "sequence_index.json"), "w"))

    man = {
        "name": "familyfm_xfam_selectivity",
        "layout": "SLS, cross family and same family" if by_type else "SLS, cross family",
        "question": ("one ligand, two targets: which target the compound prefers. "
                     "Cross family (no shared family label) and same family (two "
                     "different targets sharing one) are reported separately"
                     if by_type else "one ligand, two targets from different "
                     "families: which target the compound prefers"),
        "feature_order": "[ sequence A 480 | ligand 1038 | sequence B 480 ] = 1998",
        "source": "ChEMBL 37",
        "training": {"pairs_file": os.path.basename(a.pairs), "fit_comparisons": int(len(fit)),
                     "rows_after_swap": int(2 * len(fit)),
                     "label_balance": round(float(np.concatenate([y, 1 - y]).mean()), 4),
                     "fit_ligands": int(fit.inchikey.nunique()),
                     "targets": int(len(seqs)), "fit_seconds": round(fit_s)},
        "split_basis": "compound-disjoint by InChIKey, blake2b, %.0f%% held out"
                       % (100 * a.test_fraction),
        "hyperparameters": {"n_estimators": a.trees,
                            "min_samples_leaf": a.min_samples_leaf,
                            "max_features": "sqrt", "random_state": a.seed},
        "MEASURED_PERFORMANCE": pooled,
        "environment": {"python": sys.version.split()[0],
                        "numpy": np.__version__,
                        "sklearn": __import__("sklearn").__version__},
    }
    if by_type:
        man["MEASURED_PERFORMANCE"]["WARNING"] = (
            "these figures pool cross-family and same-family comparisons and must "
            "never be quoted. Use MEASURED_PERFORMANCE_BY_PAIR_TYPE.")
        man["MEASURED_PERFORMANCE_BY_PAIR_TYPE"] = by_type
    json.dump(man, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=1)
    cols = ["inchikey", "uniprot_a", "uniprot_b", "family_a", "family_b",
            "assay_measure", "pic50_a", "pic50_b", "p", "correct"]
    if by_type:
        cols.insert(5, "pair_type")
    tst[cols].to_csv(
        os.path.join(a.out, "heldout_scored.csv"), index=False)

    for nm, blk in (by_type or {"all": pooled}).items():
        print("\nstrength bands, %s:" % nm)
        for k, v in blk["strength_cumulative"].items():
            print("  >= %s  coverage %.3f  accuracy %s  n %d"
                  % (k, v["coverage"], v["accuracy"], v["n"]))
    print("\ngap strata:", json.dumps(gap_strata))
    print("\ntop family pairs by held-out size:")
    for k, v in sorted(per_pair.items(), key=lambda x: -x[1]["n"])[:12]:
        print("  %-60s n %6d  acc %.3f" % (k[:58], v["n"], v["accuracy"]))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
