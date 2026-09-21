"""The two baselines for the target-preference model, measured per pair type.

A held-out accuracy means little without what a simpler model reaches on the
same comparisons. Two baselines, on the same compound-disjoint split and with the
same forest settings as train.py:

  FAMILY PRIOR   for each family pairing, always pick whichever family won more
                 often in the fit set. Inside one family this is chance by
                 construction, so every point above 0.50 there is earned.
  LIGAND BLIND   the same forest fitted on the two sequence blocks alone.

Both are reported separately for cross-family and same-family comparisons,
because a pooled baseline over both arms describes neither.

    python -m familyfm.baselines_by_pair_type \\
        --pairs data/pairs/all_both.csv --bundle models/xfam_v2_samefamily \\
        --out baselines_by_pair_type.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import embed as emb                                            # noqa: E402
from train import test_side                                    # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="pairs file with a pair_type column")
    ap.add_argument("--bundle", required=True, help="the fitted bundle, for its vectors")
    ap.add_argument("--out", required=True)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--trees", type=int, default=300)
    ap.add_argument("--min-samples-leaf", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-family-n", type=int, default=200)
    a = ap.parse_args(argv)

    d = pd.read_csv(a.pairs, usecols=["sequence_a", "sequence_b", "pic50_a", "pic50_b",
                                      "family_a", "family_b", "pair_type", "inchikey"])
    d["is_test"] = [test_side(k, a.test_fraction) for k in d.inchikey]
    fit, tst = d[~d.is_test], d[d.is_test]

    idx = json.load(open(os.path.join(a.bundle, "sequence_index.json")))["seq_key_to_row"]
    V = np.load(os.path.join(a.bundle, "sequence_vectors.npz"))["vectors"]

    def sv(col):
        return V[[idx[emb.sequence_key(s)] for s in col]]

    y = (fit.pic50_a.values > fit.pic50_b.values).astype(int)
    lab = tst.pic50_a.values > tst.pic50_b.values

    # family prior: the family that won more often in the fit set
    def fk(x, z):
        return tuple(sorted((x.split("|")[0], z.split("|")[0])))
    win, tot = collections.Counter(), collections.Counter()
    for x, z, ya in zip(fit.family_a, fit.family_b, y):
        k = fk(x, z)
        tot[k] += 1
        if (ya == 1) == (x.split("|")[0] == k[0]):
            win[k] += 1
    pred = []
    for x, z in zip(tst.family_a, tst.family_b):
        k = fk(x, z)
        p = win[k] / tot[k] if tot[k] else 0.5
        pred.append((p >= 0.5) == (x.split("|")[0] == k[0]))
    prior = np.array(pred) == lab

    # ligand blind, both orders averaged exactly as the model is
    fa, fb, ta, tb = sv(fit.sequence_a), sv(fit.sequence_b), sv(tst.sequence_a), sv(tst.sequence_b)
    rb = RandomForestClassifier(n_estimators=a.trees, min_samples_leaf=a.min_samples_leaf,
                                max_features="sqrt", random_state=a.seed, n_jobs=-1)
    rb.fit(np.vstack([np.hstack([fa, fb]), np.hstack([fb, fa])]), np.concatenate([y, 1 - y]))
    cb = int(np.where(rb.classes_ == 1)[0][0])
    pb = 0.5 * (rb.predict_proba(np.hstack([ta, tb]))[:, cb]
                + 1 - rb.predict_proba(np.hstack([tb, ta]))[:, cb])
    blind = (pb > 0.5) == lab

    t = tst.pair_type.values
    out = {}
    for nm in ("cross", "same"):
        m = t == nm
        out[nm] = {"n": int(m.sum()), "family_prior": round(float(prior[m].mean()), 4),
                   "ligand_blind": round(float(blind[m].mean()), 4)}
    # same-family ligand blind by the family the two targets share
    shared = [set(x.split("|")) & set(z.split("|")) for x, z in zip(tst.family_a, tst.family_b)]
    fams = collections.Counter(f for s, k in zip(shared, t) if k == "same" for f in s)
    out["same_by_family"] = {}
    for f, n in fams.most_common():
        if n < a.min_family_n:
            continue
        m = np.array([k == "same" and f in s for s, k in zip(shared, t)])
        out["same_by_family"][f] = {"n": int(m.sum()),
                                    "ligand_blind": round(float(blind[m].mean()), 4)}
    print(json.dumps(out, indent=1))
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
