"""Inference API shipped inside every FamilyFM bundle.

One question: given one ligand and two protein targets from different families,
which target the compound prefers. The answer is a probability, never an affinity.

    import importlib.util, sys
    B = "ffm-models/xfam_v1"        # wherever install.sh put the bundle
    spec = importlib.util.spec_from_file_location("fp", B + "/predict.py")
    fp = importlib.util.module_from_spec(spec); sys.modules["fp"] = fp
    sys.path.insert(0, B); spec.loader.exec_module(fp)
    m = fp.load(B)
    fp.compare_targets(m, smiles, "P00533", "P08684")      # EGFR against CYP3A4
    fp.rank_targets(m, smiles, ["P00533", "P08684", "P31645"])

Both target orders are scored and averaged, so compare_targets(a,b) +
compare_targets(b,a) is exactly 1 and a target compared with itself returns
0.500000. Always read the band with the number: `strength` and the accuracy the
manifest records for that band.
"""
import json
import os
import sys

import numpy as np

LIG_DIMS = 1038
SEQ_DIMS = 480
SLS_DIMS = 2 * SEQ_DIMS + LIG_DIMS
_ELEMENTS = ("C", "N", "O", "S", "F", "Cl", "Br", "I", "P")


def load(bundle_dir="."):
    import joblib
    man = json.load(open(os.path.join(bundle_dir, "MANIFEST.json")))
    idx = json.load(open(os.path.join(bundle_dir, "sequence_index.json")))
    vec = np.load(os.path.join(bundle_dir, "sequence_vectors.npz"))["vectors"]
    acc = json.load(open(os.path.join(bundle_dir, "target_index.json")))
    rf = joblib.load(os.path.join(bundle_dir, "familyfm_selectivity.joblib"))
    if rf.n_features_in_ != SLS_DIMS:
        raise RuntimeError("manifest says SLS (%d features) but the forest has %d"
                           % (SLS_DIMS, rf.n_features_in_))
    return {"manifest": man, "rf": rf, "vectors": vec, "dir": bundle_dir,
            "seq_key_to_row": idx["seq_key_to_row"],
            "accession_to_seq_key": acc["accession_to_seq_key"],
            "family_of": acc.get("family_of", {}),
            "c1": int(np.where(rf.classes_ == 1)[0][0])}


def targets(m):
    """Every accession this bundle can resolve."""
    return sorted(m["accession_to_seq_key"])


def family_of(m, accession):
    return m["family_of"].get(accession)


def ligand_features(smiles_list):
    """Morgan COUNT fingerprint r=2/1024, then the 14 descriptors, in order."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    bad = [s for s in smiles_list if Chem.MolFromSmiles(s) is None]
    if bad:
        raise ValueError("RDKit could not parse %d SMILES; they are not silently "
                         "dropped. First few: %s" % (len(bad), "; ".join(bad[:5])))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    rows = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        fp = np.asarray(gen.GetCountFingerprintAsNumPy(mol), dtype=np.float32)
        c = {e: 0 for e in _ELEMENTS}
        for at in mol.GetAtoms():
            if at.GetSymbol() in c:
                c[at.GetSymbol()] += 1
        rows.append(np.concatenate([fp, np.array(
            [Descriptors.MolWt(mol), mol.GetNumHeavyAtoms(), mol.GetNumBonds(),
             Descriptors.NumRotatableBonds(mol), Descriptors.RingCount(mol)]
            + [c[e] for e in _ELEMENTS], dtype=np.float32)]))
    return np.stack(rows)


def seq_vector(m, accession):
    key = m["accession_to_seq_key"].get(accession) or \
        m["accession_to_seq_key"].get(str(accession).upper())
    if key is None:
        raise KeyError("%r is not in this bundle. predict.targets() lists what is."
                       % accession)
    return m["vectors"][m["seq_key_to_row"][key]]


def _rows(first, middle, last):
    X = np.empty((len(middle), SLS_DIMS), dtype=np.float32)
    X[:, :SEQ_DIMS] = first
    X[:, SEQ_DIMS:SEQ_DIMS + LIG_DIMS] = middle
    X[:, SEQ_DIMS + LIG_DIMS:] = last
    return X


def _both_orders(m, first, middle, last):
    """THE aggregation: score each comparison both ways and average."""
    p1 = m["rf"].predict_proba(_rows(first, middle, last))[:, m["c1"]]
    p2 = m["rf"].predict_proba(_rows(last, middle, first))[:, m["c1"]]
    return 0.5 * (p1 + (1.0 - p2))


def strength(p):
    """The larger of the two probabilities. Read the manifest for its accuracy."""
    return float(max(p, 1.0 - p))


def compare_targets(m, smiles, accession_a, accession_b):
    """-> P(this ligand prefers target A over target B)."""
    F = ligand_features([smiles])
    va = np.asarray(seq_vector(m, accession_a), dtype=np.float32)[None, :]
    vb = np.asarray(seq_vector(m, accession_b), dtype=np.float32)[None, :]
    return float(_both_orders(m, va, F, vb)[0])


def compare_many_targets(m, smiles_list, accession_a, accession_b):
    """Many ligands, one pair of targets, in a single forest call."""
    smiles_list = list(smiles_list)
    if not smiles_list:
        return np.zeros(0)
    F = ligand_features(smiles_list)
    n = len(smiles_list)
    va = np.repeat(np.asarray(seq_vector(m, accession_a),
                              dtype=np.float32)[None, :], n, axis=0)
    vb = np.repeat(np.asarray(seq_vector(m, accession_b),
                              dtype=np.float32)[None, :], n, axis=0)
    return _both_orders(m, va, F, vb)


def rank_targets(m, smiles, accessions):
    """-> [(accession, family, mean win rate, mean strength)], best first.

    Every target is compared with every other, both orders averaged, and the win
    rate is the mean over those comparisons. With one ligand and n targets this is
    n(n-1)/2 comparisons.
    """
    accs = list(dict.fromkeys(accessions))
    if len(accs) < 2:
        raise ValueError("rank_targets needs at least two targets")
    F = ligand_features([smiles])
    V = {a: np.asarray(seq_vector(m, a), dtype=np.float32) for a in accs}
    win = {a: [] for a in accs}
    st = {a: [] for a in accs}
    for i in range(len(accs)):
        for j in range(i + 1, len(accs)):
            a, b = accs[i], accs[j]
            p = float(_both_orders(m, V[a][None, :], F, V[b][None, :])[0])
            win[a].append(p)
            win[b].append(1.0 - p)
            st[a].append(strength(p))
            st[b].append(strength(p))
    out = [(a, family_of(m, a), float(np.mean(win[a])), float(np.mean(st[a])))
           for a in accs]
    return sorted(out, key=lambda r: -r[2])
