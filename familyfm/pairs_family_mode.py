"""Measurements to target-preference comparisons, across and within families.

One ligand, two targets, and the question is which target the compound prefers.
Every row is stamped with `pair_type`: `cross` when the two targets share no
family label, `same` when two different targets share one. The two kinds are
fitted together and REPORTED SEPARATELY (see train.py); a single pooled accuracy
over both would be the one number nobody can act on.

The rules:

  ENDPOINT MATCHED   both readings carry the same endpoint, so a Ki is never
                     ranked against an IC50
  TIES DROPPED       equal readings carry no answer
  SIDE RANDOMIZED    which target is A is a coin flip per comparison
  CAPS               one per ligand and pair type, one per (family group,
                     endpoint)
  ONE PROTEIN        a pair whose two entries are the same UniProt accession is
                     never formed. 45 human proteins in ChEMBL carry two level-2
                     labels, and that is the case the rule exists for

`--mode` selects which kinds are formed:

  both    everything, each row labeled (the released model)
  cross   only targets that share no family label (what pairs.py forms)
  same    only different targets that share a family label

    python -m familyfm.pairs_family_mode \
        --measurements data/measurements/measurements_all.csv \
        --out data/pairs/all_both.csv --mode both
"""
import argparse
import collections
import csv
import itertools
import json
import os

import numpy as np
import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--measurements", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=("cross", "same", "both"), default="both")
    ap.add_argument("--per-ligand-cap", type=int, default=50)
    ap.add_argument("--per-pair-cap", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    d = pd.read_csv(a.measurements, low_memory=False)
    print("measurements: %d" % len(d), flush=True)
    seq_of = dict(zip(d.uniprot, d.sequence))
    fam_of = {u: frozenset(str(f).split("|")) for u, f in zip(d.uniprot, d.family)}
    smi_of = dict(zip(d.inchikey, d.smiles))

    groups = collections.defaultdict(dict)
    for k, u, p, e in zip(d.inchikey, d.uniprot, d.pic50, d.assay_measure):
        groups[(k, e)][u] = p

    cand = collections.defaultdict(list)
    per_lig = collections.Counter()
    formed = collections.Counter()
    ties = collections.Counter()
    for (key, end), tg in groups.items():
        if len(tg) < 2:
            continue
        us = list(tg)
        pool = []
        for i, j in itertools.combinations(range(len(us)), 2):
            ua, ub = us[i], us[j]
            if ua == ub:                      # never a comparison
                continue
            shared = fam_of[ua] & fam_of[ub]
            kind = "same" if shared else "cross"
            if a.mode != "both" and kind != a.mode:
                continue
            if tg[ua] == tg[ub]:
                ties[kind] += 1
                continue
            pool.append((ua, ub, kind))
        if not pool:
            continue
        # THE CAP IS PER PAIR TYPE, NOT SHARED. Same-family candidates outnumber
        # cross-family ones about 3 to 1, and a shared 50-slot budget let them
        # crowd the cross-family comparisons down from 89,888 to 73,871. Each
        # kind gets its own budget, so the cross-family arm keeps the SIZE that
        # pairs.py produces. Its membership is not identical: where a compound
        # carries more than 50 candidates of one kind, which 50 are drawn comes
        # from one random stream that in --mode both also serves the same-family
        # draws, so it diverges from pairs.py. Measured: 8,446 of 89,888
        # cross-family pairs differ.
        kept = []
        for kind in ("cross", "same"):
            sub = [q for q in pool if q[2] == kind]
            if len(sub) > a.per_ligand_cap:
                pick = rng.choice(len(sub), a.per_ligand_cap, replace=False)
                sub = [sub[t] for t in pick]
            kept.extend(sub)
        pool = kept
        for ua, ub, kind in pool:
            # ONE ROW PER PAIR, side chosen by a coin flip. This is deliberate
            # and must not be "improved" by emitting both orientations:
            # train.py already does swap augmentation at fit time, stacking
            # (A,B) with (B,A) and inverting the label, which is why its label
            # balance is exactly 0.5000 and compare(A,B)+compare(B,A)==1. A
            # pre-doubled file would be augmented again and every pair would
            # enter four times.
            if rng.random() < 0.5:
                ua, ub = ub, ua
            fa, fb = sorted(fam_of[ua])[0], sorted(fam_of[ub])[0]
            cand[(tuple(sorted((fa, fb))), end, kind)].append(
                {"smiles": smi_of[key], "sequence_a": seq_of[ua],
                 "sequence_b": seq_of[ub], "pic50_a": tg[ua], "pic50_b": tg[ub],
                 "relation_a": "=", "relation_b": "=", "assay_measure": end,
                 "uniprot_a": ua, "uniprot_b": ub,
                 "family_a": "|".join(sorted(fam_of[ua])),
                 "family_b": "|".join(sorted(fam_of[ub])),
                 "pair_type": kind, "inchikey": key})
            formed[kind] += 1
            per_lig[key] += 1

    rows, capped = [], 0
    for k, v in cand.items():
        if len(v) > a.per_pair_cap:
            pick = rng.choice(len(v), a.per_pair_cap, replace=False)
            v = [v[t] for t in pick]
            capped += 1
        rows.extend(v)
    rng.shuffle(rows)
    cols = ["smiles", "sequence_a", "sequence_b", "pic50_a", "pic50_b",
            "relation_a", "relation_b", "assay_measure", "uniprot_a",
            "uniprot_b", "family_a", "family_b", "pair_type", "inchikey"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    by_type = collections.Counter(r["pair_type"] for r in rows)
    meta = {"mode": a.mode, "comparisons": len(rows),
            "by_pair_type": dict(by_type),
            "formed_before_caps": dict(formed), "ties_dropped": dict(ties),
            "ligands": len({r["inchikey"] for r in rows}),
            "targets": len({r["uniprot_a"] for r in rows} |
                           {r["uniprot_b"] for r in rows}),
            "groups_capped": capped,
            "per_ligand_cap": a.per_ligand_cap,
            "per_pair_cap": a.per_pair_cap,
            "A_wins_fraction": round(
                sum(1 for r in rows if r["pic50_a"] > r["pic50_b"])
                / max(len(rows), 1), 4),
            "by_endpoint": dict(collections.Counter(
                r["assay_measure"] for r in rows).most_common()),
            "by_endpoint_and_type": {
                "%s/%s" % (k[0], k[1]): v for k, v in
                collections.Counter((r["assay_measure"], r["pair_type"])
                                    for r in rows).items()}}
    json.dump(meta, open(a.out.replace(".csv", "_meta.json"), "w"), indent=1)
    print("\nwrote %s" % a.out)
    print("  comparisons %d  (%s)" % (len(rows), dict(by_type)))
    print("  ligands %d, targets %d, capped groups %d"
          % (meta["ligands"], meta["targets"], capped))
    print("  A wins %.4f (0.5 means the side is randomized)"
          % meta["A_wins_fraction"])
    print("  ties dropped %s" % dict(ties))


if __name__ == "__main__":
    main()
