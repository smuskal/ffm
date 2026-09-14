"""Fit and evaluate the cross-family target-preference comparator, with its baselines.

Split is compound-disjoint: a ligand's InChIKey decides which side it lands on, so
no ligand appears in both fit and test. Ties are already dropped at pairing. Every
comparison is entered twice with the label inverted, and every prediction averages
both target orders, so compare(A,B) + compare(B,A) is 1 and label balance is
exactly 0.5000.

Three numbers are reported together, and the first is meaningless without the
other two:

  the model
  the FAMILY PRIOR: for each family pair, always pick whichever family won more
    often in the fit set. If the model cannot beat this, it is a lookup table of
    family pairs and the ligand is doing nothing.
  LIGAND BLIND: the same forest trained on the two sequence blocks alone. This is
    the same test in model form.

    python -m familyfm.train --pairs data/pairs/all.csv --out models/xfam_v1
"""
import argparse
import collections
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--trees", type=int, default=300)
    ap.add_argument("--min-samples-leaf", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip-ligand-blind", action="store_true")
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
    print("\nMODEL held-out accuracy %.4f on %d comparisons" % (acc, len(tst)), flush=True)

    # ---- baseline 1: family prior from the fit set ----
    def famkey(r):
        return tuple(sorted((r[0].split("|")[0], r[1].split("|")[0])))
    win = collections.Counter()
    tot = collections.Counter()
    for faa, fbb, ya in zip(fit.family_a, fit.family_b, y):
        k = famkey((faa, fbb))
        tot[k] += 1
        first = faa.split("|")[0]
        # did the family that sorts first win this comparison
        if (ya == 1 and first == k[0]) or (ya == 0 and first != k[0]):
            win[k] += 1
    prior_pred = []
    for faa, fbb in zip(tst.family_a, tst.family_b):
        k = famkey((faa, fbb))
        p = win[k] / tot[k] if tot[k] else 0.5
        first_is_k0 = faa.split("|")[0] == k[0]
        # predict A wins if the family on side A is the one that usually wins
        prior_pred.append((p >= 0.5) == first_is_k0)
    prior_acc = float((np.array(prior_pred) == lab).mean())
    print("FAMILY PRIOR accuracy   %.4f" % prior_acc, flush=True)

    # ---- baseline 2: ligand blind, same forest without the ligand block ----
    blind_acc = None
    if not a.skip_ligand_blind:
        Xb = np.vstack([np.hstack([fa, fb]), np.hstack([fb, fa])])
        Yb = np.concatenate([y, 1 - y])
        rb = RandomForestClassifier(n_estimators=a.trees,
                                    min_samples_leaf=a.min_samples_leaf,
                                    max_features="sqrt", random_state=a.seed,
                                    n_jobs=-1)
        rb.fit(Xb, Yb)
        cb = int(np.where(rb.classes_ == 1)[0][0])
        pb1 = rb.predict_proba(np.hstack([ta, tb]))[:, cb]
        pb2 = rb.predict_proba(np.hstack([tb, ta]))[:, cb]
        Pb = 0.5 * (pb1 + (1.0 - pb2))
        blind_acc = float((((Pb > 0.5) == lab)).mean())
        print("LIGAND BLIND accuracy   %.4f" % blind_acc, flush=True)
        del Xb, Yb

    # ---- per family pair and per endpoint ----
    tst = tst.assign(p=P, correct=correct, strength=strength)
    per_pair = {}
    for k, g in tst.groupby([tst.family_a.str.split("|").str[0],
                             tst.family_b.str.split("|").str[0]]):
        kk = " | ".join(sorted(k))
        if kk in per_pair or len(g) < 30:
            continue
        per_pair[kk] = {"n": int(len(g)), "accuracy": round(float(g.correct.mean()), 4)}
    per_end = {k: {"n": int(len(g)), "accuracy": round(float(g.correct.mean()), 4)}
               for k, g in tst.groupby("assay_measure")}
    gap = np.abs(tst.pic50_a - tst.pic50_b)
    gap_strata = {}
    for lo, hi, nm in [(0, .5, "under 0.5 log"), (.5, 1, "0.5 to 1 log"),
                       (1, 2, "1 to 2 logs"), (2, 99, "over 2 logs")]:
        m = (gap >= lo) & (gap < hi)
        if m.sum():
            gap_strata[nm] = {"n": int(m.sum()),
                              "accuracy": round(float(tst.correct[m].mean()), 4)}

    joblib.dump(rf, os.path.join(a.out, "familyfm_selectivity.joblib"), compress=3)
    idx = {emb.sequence_key(s): i for i, s in enumerate(seqs)}
    np.savez_compressed(os.path.join(a.out, "sequence_vectors.npz"),
                        vectors=np.stack([SV[s] for s in seqs]))
    json.dump({"seq_key_to_row": idx}, open(os.path.join(a.out, "sequence_index.json"), "w"))

    man = {
        "name": "familyfm_xfam_selectivity", "layout": "SLS, cross family",
        "question": "one ligand, two targets from different families: "
                    "which target the compound prefers",
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
        "MEASURED_PERFORMANCE": {
            "held_out_comparisons": int(len(tst)),
            "held_out_ligands": int(tst.inchikey.nunique()),
            "accuracy": round(acc, 4),
            "family_prior_accuracy": round(prior_acc, 4),
            "ligand_blind_accuracy": None if blind_acc is None else round(blind_acc, 4),
            "strength_cumulative": cumulative(strength, correct),
            "by_endpoint": per_end, "by_family_pair": per_pair,
            "gap_strata": gap_strata},
        "environment": {"python": sys.version.split()[0],
                        "numpy": np.__version__,
                        "sklearn": __import__("sklearn").__version__},
    }
    json.dump(man, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=1)
    tst[["inchikey", "uniprot_a", "uniprot_b", "family_a", "family_b",
         "assay_measure", "pic50_a", "pic50_b", "p", "correct"]].to_csv(
        os.path.join(a.out, "heldout_scored.csv"), index=False)

    print("\nstrength bands:")
    for k, v in man["MEASURED_PERFORMANCE"]["strength_cumulative"].items():
        print("  >= %s  coverage %.3f  accuracy %s  n %d"
              % (k, v["coverage"], v["accuracy"], v["n"]))
    print("\ngap strata:", json.dumps(gap_strata))
    print("\ntop family pairs by held-out size:")
    for k, v in sorted(per_pair.items(), key=lambda x: -x[1]["n"])[:12]:
        print("  %-60s n %6d  acc %.3f" % (k[:58], v["n"], v["accuracy"]))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
