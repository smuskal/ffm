"""ChEMBL 37 to one measurement per row, scoped for target preference.

Built from ChEMBL 37.

Scope:

  human, target_type SINGLE PROTEIN, carrying a level-2 protein class (the family)
  no PubChem BioAssay rows (src_id 7), which are high-throughput screening
  endpoints Ki, IC50, Kd, EC50, Kb only
  pActivity derived from standard_value and standard_units, NEVER from
    pchembl_value, which in this data is present on erroneous rows
  pActivity kept in (-1, 13]
  exact readings only: standard_relation '=' or null. Censored readings are
    dropped rather than ordered, which is the strict choice for a first build
  one value per (uniprot, inchikey, endpoint) by MEDIAN, never by most potent

Output: data/measurements/measurements_all.csv with the measurement schema
plus a `family` column.

    python -m familyfm.extract --out data/measurements
"""
import argparse
import collections
import csv
import json
import os
import sqlite3
import statistics
import sys

# Path to a local ChEMBL 37 SQLite. Override with the CHEMBL_DB environment
# variable or --db; the default assumes the release is beside this checkout.
DB = os.environ.get("CHEMBL_DB", "data/chembl_37.db")
UNIT = {"M": 0, "mM": -3, "uM": -6, "nM": -9, "pM": -12, "fM": -15}
ENDPOINTS = ("Ki", "IC50", "Kd", "EC50", "Kb")
PACT_LO, PACT_HI = -1.0, 13.0


def family_map(c):
    """-> {protein_class_id: level-2 family name}."""
    rows = c.execute("select protein_class_id,parent_id,pref_name,class_level "
                     "from protein_classification").fetchall()
    par = {r[0]: r[1] for r in rows}
    name = {r[0]: r[2] for r in rows}
    lvl = {r[0]: r[3] for r in rows}

    def anc2(cid):
        seen = set()
        while cid is not None and cid not in seen:
            seen.add(cid)
            if lvl.get(cid) == 2:
                return cid
            cid = par.get(cid)
        return None
    return {cid: name.get(anc2(cid)) for cid in par}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=DB)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    c = sqlite3.connect("file:%s?mode=ro" % a.db, uri=True)

    famof = family_map(c)
    tinfo = {}                      # tid -> (accession, sequence, frozenset(families))
    fams = collections.defaultdict(set)
    seqs = {}
    accs = {}
    for tid, acc, seq, cid in c.execute(
            """select td.tid, cs.accession, cs.sequence, cc.protein_class_id
               from target_dictionary td
               join target_components tc on td.tid=tc.tid
               join component_sequences cs on tc.component_id=cs.component_id
               join component_class cc on tc.component_id=cc.component_id
               where td.organism='Homo sapiens' and td.target_type='SINGLE PROTEIN'
                 and cs.sequence is not null"""):
        f = famof.get(cid)
        if f:
            fams[tid].add(f)
            seqs[tid] = seq
            accs[tid] = acc
    for tid in fams:
        tinfo[tid] = (accs[tid], seqs[tid], frozenset(fams[tid]))
    print("targets in scope: %d" % len(tinfo), flush=True)

    # ---- activities ----
    vals = collections.defaultdict(list)     # (acc, molregno, endpoint) -> [pAct]
    tgt_of = {}
    drop = collections.Counter()
    n = 0
    q = """select a.molregno, ass.tid, a.standard_type, a.standard_value,
                  a.standard_units, a.standard_relation
           from activities a join assays ass on a.assay_id=ass.assay_id
           where ass.src_id != 7 and a.standard_value is not null
             and a.standard_type in ('Ki','IC50','Kd','EC50','Kb')"""
    for mol, tid, st, val, unit, rel in c.execute(q):
        n += 1
        info = tinfo.get(tid)
        if info is None:
            drop["target out of scope"] += 1
            continue
        if unit not in UNIT:
            drop["unit not a concentration"] += 1
            continue
        if rel not in (None, "", "="):
            drop["censored reading"] += 1
            continue
        try:
            v = float(val)
        except (TypeError, ValueError):
            drop["value not numeric"] += 1
            continue
        if v <= 0:
            drop["value not positive"] += 1
            continue
        import math
        p = -math.log10(v * (10.0 ** UNIT[unit]))
        if not (PACT_LO < p <= PACT_HI):
            drop["pActivity out of range"] += 1
            continue
        acc = info[0]
        vals[(acc, mol, st)].append(p)
        tgt_of[acc] = tid
        if n % 1000000 == 0:
            print("  %d rows read" % n, flush=True)
    print("rows read: %d" % n, flush=True)
    for k, v in drop.most_common():
        print("  dropped, %-28s %d" % (k, v), flush=True)

    # ---- structures for the molecules we kept ----
    keep = {k[1] for k in vals}
    print("molecules to resolve: %d" % len(keep), flush=True)
    smi = {}
    key = {}
    for mol, s, ik in c.execute("select molregno, canonical_smiles, "
                                "standard_inchi_key from compound_structures"):
        if mol in keep and s and ik:
            smi[mol] = s
            key[mol] = ik
    print("molecules with a structure: %d" % len(smi), flush=True)

    # ---- collapse and write ----
    path = os.path.join(a.out, "measurements_all.csv")
    rows = 0
    per_fam = collections.Counter()
    per_end = collections.Counter()
    tgts = set()
    ligs = set()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["smiles", "sequence", "pic50", "relation", "uniprot",
                    "family", "assay_measure", "inchikey"])
        for (acc, mol, st), v in vals.items():
            if mol not in smi:
                continue
            tid = tgt_of[acc]
            _, seq, ff = tinfo[tid]
            med = statistics.median(v)
            w.writerow([smi[mol], seq, round(med, 4), "=", acc,
                        "|".join(sorted(ff)), st, key[mol]])
            rows += 1
            for f in ff:
                per_fam[f] += 1
            per_end[st] += 1
            tgts.add(acc)
            ligs.add(key[mol])
    meta = {"source": "ChEMBL 37, %s" % a.db,
            "scope": "human SINGLE PROTEIN with a level-2 protein class",
            "endpoints": list(ENDPOINTS),
            "screening_excluded": "assays.src_id = 7 (PubChem BioAssay)",
            "censored": "dropped; exact readings only",
            "pactivity": "derived from standard_value and standard_units; "
                         "pchembl_value not consulted",
            "aggregation": "median per (uniprot, inchikey, endpoint)",
            "measurements": rows, "targets": len(tgts), "ligands": len(ligs),
            "by_endpoint": dict(per_end.most_common()),
            "by_family": dict(per_fam.most_common()),
            "rows_read": n, "dropped": dict(drop)}
    json.dump(meta, open(os.path.join(a.out, "measurements_meta.json"), "w"),
              indent=1)
    print("\nwrote %s" % path)
    print("  measurements %d over %d targets and %d ligands"
          % (rows, len(tgts), len(ligs)))
    print("  by endpoint: %s" % dict(per_end.most_common()))
    print("  top families: %s" % dict(per_fam.most_common(8)))


if __name__ == "__main__":
    main()
