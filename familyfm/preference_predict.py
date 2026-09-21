"""Inference API for the preference bundle.

One question: given one target and two compounds, which compound the target
prefers, that is which is active at the lower concentration there. The answer is
a probability, never a value on the potency scale.

    import familyfm.preference_predict as P
    m = P.load("ffm-models/lsl_v2_stratified")
    P.compare_ligands(m, smiles_a, smiles_b, "P00533")     # at EGFR
    P.rank_ligands(m, [s1, s2, s3, s4], "P00533")          # best first

The identical file ships inside the bundle as predict.py, so a bundle can be
loaded by file path with nothing installed.

Both compound orders are scored and averaged, so compare_ligands(a,b) +
compare_ligands(b,a) is exactly 1 and a compound against itself returns
0.500000. Always read the band with the number: `strength` and the accuracy the
manifest records for that band.
"""
import json
import os
import sys

import numpy as np

LIG_DIMS = 1038
SEQ_DIMS = 480
LSL_DIMS = 2 * LIG_DIMS + SEQ_DIMS
_ELEMENTS = ("C", "N", "O", "S", "F", "Cl", "Br", "I", "P")


def load(bundle_dir="."):
    import joblib
    man = json.load(open(os.path.join(bundle_dir, "MANIFEST.json")))
    idx = json.load(open(os.path.join(bundle_dir, "sequence_index.json")))
    vec = np.load(os.path.join(bundle_dir, "sequence_vectors.npz"))["vectors"]
    acc = json.load(open(os.path.join(bundle_dir, "target_index.json")))
    rf = joblib.load(os.path.join(bundle_dir, "familyfm_preference.joblib"))
    if rf.n_features_in_ != LSL_DIMS:
        raise RuntimeError("manifest says LSL (%d features) but the forest has %d"
                           % (LSL_DIMS, rf.n_features_in_))
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
    X = np.empty((len(middle), LSL_DIMS), dtype=np.float32)
    X[:, :LIG_DIMS] = first
    X[:, LIG_DIMS:LIG_DIMS + SEQ_DIMS] = middle
    X[:, LIG_DIMS + SEQ_DIMS:] = last
    return X


def _both_orders(m, first, middle, last):
    """THE aggregation: score each comparison both ways and average."""
    p1 = m["rf"].predict_proba(_rows(first, middle, last))[:, m["c1"]]
    p2 = m["rf"].predict_proba(_rows(last, middle, first))[:, m["c1"]]
    return 0.5 * (p1 + (1.0 - p2))


def strength(p):
    """The larger of the two probabilities. Read the manifest for its accuracy."""
    return float(max(p, 1.0 - p))


def compare_ligands(m, smiles_a, smiles_b, accession):
    """-> P(the target prefers compound A over compound B)."""
    F = ligand_features([smiles_a, smiles_b])
    v = np.asarray(seq_vector(m, accession), dtype=np.float32)[None, :]
    return float(_both_orders(m, F[0:1], v, F[1:2])[0])


def compare_many_ligands(m, smiles_list, reference, accession):
    """Many candidates against one reference compound, in a single forest call."""
    smiles_list = list(smiles_list)
    if not smiles_list:
        return np.zeros(0)
    F = ligand_features(smiles_list + [reference])
    n = len(smiles_list)
    v = np.repeat(np.asarray(seq_vector(m, accession),
                             dtype=np.float32)[None, :], n, axis=0)
    return _both_orders(m, F[:-1], v, np.repeat(F[-1:], n, axis=0))


def rank_ligands(m, smiles_list, accession):
    """-> [(smiles, mean win rate, mean strength)], the target's preference order.

    Every compound is compared with every other, both orders averaged, so n
    compounds cost n(n-1)/2 comparisons.
    """
    smis = list(dict.fromkeys(smiles_list))
    if len(smis) < 2:
        raise ValueError("rank_ligands needs at least two compounds")
    F = ligand_features(smis)
    v = np.asarray(seq_vector(m, accession), dtype=np.float32)
    pairs = [(i, j) for i in range(len(smis)) for j in range(i + 1, len(smis))]
    P = _both_orders(m, F[[i for i, _ in pairs]],
                     np.repeat(v[None, :], len(pairs), axis=0),
                     F[[j for _, j in pairs]])
    win = {s: [] for s in smis}
    st = {s: [] for s in smis}
    for (i, j), p in zip(pairs, P):
        win[smis[i]].append(float(p))
        win[smis[j]].append(1.0 - float(p))
        st[smis[i]].append(strength(p))
        st[smis[j]].append(strength(p))
    return sorted([(s, float(np.mean(win[s])), float(np.mean(st[s]))) for s in smis],
                  key=lambda r: -r[1])
