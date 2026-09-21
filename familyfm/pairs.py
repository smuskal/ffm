"""Measurements to cross-family SLS comparisons.

One ligand, two targets, and the question is which target the compound prefers.
Four rules make the comparison fair, and all four came out of the kinase and GPCR
builds:

  ENDPOINT MATCHED. Both readings carry the same endpoint, so a Ki is never
  ranked against an IC50.

  DISJOINT FAMILIES. The two targets share no family label. 45 human proteins in
  ChEMBL carry two family labels, reader and writer most often, and without this
  rule one protein filed twice would pose as a cross-family pair.

  TIES DROPPED. Equal readings carry no answer, and entering one anyway teaches
  the forest that whichever side is written first wins.

  SIDE RANDOMIZED. Which target is A is decided by a coin flip per comparison, and
  the A-wins fraction is reported. If it drifts from 0.5 the forest can score well
  by reading column position rather than chemistry.

Two caps stop a handful of heavily profiled compounds becoming the model: one per
ligand, one per family pair and endpoint.

    python -m familyfm.pairs --measurements data/measurements/measurements_all.csv \
        --out data/pairs/all.csv
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
    ap.add_argument("--per-ligand-cap", type=int, default=50)
    ap.add_argument("--per-pair-cap", type=int, default=20000,
                    help="cap per (family pair, endpoint)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    d = pd.read_csv(a.measurements, low_memory=False)
    print("measurements: %d" % len(d), flush=True)
    seq_of = dict(zip(d.uniprot, d.sequence))
    fam_of = {u: frozenset(f.split("|")) for u, f in zip(d.uniprot, d.family)}
    smi_of = dict(zip(d.inchikey, d.smiles))

    groups = collections.defaultdict(dict)          # (inchikey, endpoint) -> {uniprot: pic50}
    for k, u, p, e in zip(d.inchikey, d.uniprot, d.pic50, d.assay_measure):
        groups[(k, e)][u] = p

    cand = collections.defaultdict(list)            # (fam pair, endpoint) -> rows
    per_lig = collections.Counter()
    formed = ties = 0
    for (key, end), tg in groups.items():
        if len(tg) < 2:
            continue
        us = list(tg)
        pool = []
        for i, j in itertools.combinations(range(len(us)), 2):
            ua, ub = us[i], us[j]
            if fam_of[ua] & fam_of[ub]:
                continue
            if tg[ua] == tg[ub]:
                ties += 1
                continue
            pool.append((ua, ub))
        if not pool:
            continue
        if len(pool) > a.per_ligand_cap:
            pick = rng.choice(len(pool), a.per_ligand_cap, replace=False)
            pool = [pool[t] for t in pick]
        for ua, ub in pool:
            if rng.random() < 0.5:                  # randomize which side is A
                ua, ub = ub, ua
            fa, fb = sorted(fam_of[ua])[0], sorted(fam_of[ub])[0]
            cand[(tuple(sorted((fa, fb))), end)].append(
                {"smiles": smi_of[key], "sequence_a": seq_of[ua],
                 "sequence_b": seq_of[ub], "pic50_a": tg[ua], "pic50_b": tg[ub],
                 "relation_a": "=", "relation_b": "=", "assay_measure": end,
                 "uniprot_a": ua, "uniprot_b": ub,
                 "family_a": "|".join(sorted(fam_of[ua])),
                 "family_b": "|".join(sorted(fam_of[ub])), "inchikey": key})
            formed += 1
            per_lig[key] += 1

    rows = []
    capped = 0
    for k, v in cand.items():
        if len(v) > a.per_pair_cap:
            pick = rng.choice(len(v), a.per_pair_cap, replace=False)
            v = [v[t] for t in pick]
            capped += 1
        rows.extend(v)
    rng.shuffle(rows)
    cols = ["smiles", "sequence_a", "sequence_b", "pic50_a", "pic50_b",
            "relation_a", "relation_b", "assay_measure", "uniprot_a", "uniprot_b",
            "family_a", "family_b", "inchikey"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    a_wins = sum(1 for r in rows if r["pic50_a"] > r["pic50_b"]) / max(len(rows), 1)
    fp = collections.Counter()
    for r in rows:
        fp[tuple(sorted((r["family_a"].split("|")[0],
                         r["family_b"].split("|")[0])))] += 1
    meta = {"comparisons": len(rows), "formed_before_caps": formed,
            "ties_dropped": ties, "ligands": len({r["inchikey"] for r in rows}),
            "targets": len({r["uniprot_a"] for r in rows} |
                           {r["uniprot_b"] for r in rows}),
            "family_pairs": len(fp), "groups_capped": capped,
            "per_ligand_cap": a.per_ligand_cap, "per_pair_cap": a.per_pair_cap,
            "A_wins_fraction": round(a_wins, 4),
            "by_endpoint": dict(collections.Counter(
                r["assay_measure"] for r in rows).most_common()),
            "top_family_pairs": [{"a": k[0], "b": k[1], "comparisons": v}
                                 for k, v in fp.most_common(20)]}
    json.dump(meta, open(a.out.replace(".csv", "_meta.json"), "w"), indent=1)
    print("\nwrote %s" % a.out)
    print("  comparisons %d on %d ligands, %d targets, %d family pairs"
          % (len(rows), meta["ligands"], meta["targets"], meta["family_pairs"]))
    print("  ties dropped %d | A wins %.4f (0.5 means side is randomized)"
          % (ties, a_wins))
    print("  by endpoint %s" % meta["by_endpoint"])
    for k, v in fp.most_common(10):
        print("    %-32s %-32s %7d" % (k[0][:30], k[1][:30], v))


if __name__ == "__main__":
    main()
