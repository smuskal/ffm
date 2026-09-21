"""The feature recipe. Self-contained: nothing here imports anything else.

Ligand block  1,038 dims = Morgan COUNT fingerprint radius 2 at 1,024 bits
              (rdFingerprintGenerator, NOT GetMorganFingerprintAsBitVect) then 14
              descriptors in this exact order: MolWt, HeavyAtomCount, NumBonds,
              NumRotatableBonds, RingCount, then counts of C N O S F Cl Br I P.
Sequence block  480 dims = ESM2 facebook/esm2_t12_35M_UR50D, mean-pooled over
              residues with BOS and EOS excluded, first 480 dims. See embed.py.

Target-preference width, SLS = 480 + 1038 + 480 = 1,998, sequence, ligand, sequence.
Compound-preference width, LSL = 1038 + 480 + 1038 = 2,556, ligand, sequence, ligand.

`rdkit_check.py` in this directory verifies these features match the ones the
shipped forest was fitted on. Run it before trusting a build.
"""
import numpy as np

LIG_DIMS = 1038
SEQ_DIMS = 480
SLS_DIMS = 2 * SEQ_DIMS + LIG_DIMS
_ELEMENTS = ("C", "N", "O", "S", "F", "Cl", "Br", "I", "P")


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
    F = np.stack(rows)
    if F.shape[1] != LIG_DIMS:
        raise ValueError("ligand features are %d wide, expected %d"
                         % (F.shape[1], LIG_DIMS))
    return F


def sls_rows(seq_a, lig, seq_b):
    """[ sequence A 480 | ligand 1038 | sequence B 480 ], the SLS feature order."""
    n = len(lig)
    X = np.empty((n, SLS_DIMS), dtype=np.float32)
    X[:, :SEQ_DIMS] = seq_a
    X[:, SEQ_DIMS:SEQ_DIMS + LIG_DIMS] = lig
    X[:, SEQ_DIMS + LIG_DIMS:] = seq_b
    return X
