"""Prove two environments compute the SAME ligand features, or find where they differ.

RDKit computes 1,024 of every ligand block's 1,038 dimensions. A service built on
a different RDKit than the forest was fitted on produces different features, scores
worse, and nothing in the bundle notices. Asserting "we both use 2025.09.5" is not
proof; this is.

    python -m familyfm.rdkit_check --write    once, to record the reference
    python -m familyfm.rdkit_check            anywhere else, to verify

The check hashes the exact float32 matrix the model is fed, so it catches a version
difference, a different fingerprint API, a platform float difference, and a
sanitization change alike.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
# the reference ships with the repository, so verify mode works on a fresh
# clone with nothing built yet
DEST = os.path.join(PROJECT, "docs", "rdkit_fingerprint_check.json")

# Structures chosen to exercise the encoder: aromatics, stereocenters, charges,
# tautomer-prone rings, a macrocycle, a peptide, and explicit hydrogen isotopes.
SMILES = [
    "CCO",
    "c1ccccc1",
    "CC(=O)Nc1ccc(O)cc1",
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "C[C@@H](C(=O)O)N",
    "C[C@H](C(=O)O)N",
    "OC(=O)c1ccccc1O",
    "c1ccc2[nH]ccc2c1",
    "C1CCNCC1",
    "[Na+].[Cl-]",
    "CC(=O)[O-]",
    "c1cc[nH+]cc1",
    "O=S(=O)(N)c1ccc(Cl)cc1",
    "FC(F)(F)c1ccccc1",
    "Clc1ccc(Br)cc1I",
    "CCCCCCCCCCCCCCCC(=O)O",
    "N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(=O)O)C(=O)O",
    "C1CCC(CC1)N2CCN(CC2)c3ncccn3",
    ("Cc1cc(-n2nc3c(c2-n2ccn(-c4ccc5c(cnn5C)c4F)c2=O)[C@H](C)N(C(=O)"
     "c2cc4cc([C@H]5CCOC(C)(C)C5)ccc4n2[C@@]2(c4noc(=O)[nH]4)C[C@@H]2C)"
     "CC3)cc(C)c1F"),
    # Explicit hydrogen isotopes. RDKit 2025.09.6 stopped counting a fully
    # deuterated methyl as a rotatable bond, so NumRotatableBonds moves for these
    # and not for their protiated analogues. Without them in this set the check
    # passes across that boundary and the change reaches a model silently.
    "[2H]C([2H])([2H])N(C)C(=O)c1ccccc1",
    "[2H]C([2H])([2H])CO",
    ("[2H]C([2H])([2H])N(C(=O)c1c(F)cccc1Cl)c1ccc(-c2cc(NC(C)=O)nn2C(C)C)"
     "cc1N1CCCCC1"),
    "[3H]CC(=O)Nc1ccccc1",
]


def features():
    """The encoder, reached both ways a user reaches it.

    The scoring path (familyfm.features) and the published extension tool
    (extend/ffm_extend.py) each carry the recipe; they must agree with each other
    as well as with the reference, so a drift in either is caught.
    """
    from familyfm.features import ligand_features
    sys.path.insert(0, os.path.join(PROJECT, "extend"))
    from ffm_extend import ligand_matrix
    F = np.asarray(ligand_features(SMILES), dtype=np.float32)
    G = np.asarray(ligand_matrix(SMILES), dtype=np.float32)
    if F.shape != G.shape or not np.array_equal(F, G):
        sys.exit("MISMATCH. familyfm.features and extend/ffm_extend.py compute "
                 "different ligand features in this environment.")
    return F


def versions():
    import rdkit.rdBase
    import sklearn
    return {"python": sys.version.split()[0],
            "rdkit": rdkit.rdBase.rdkitVersion,
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
            "platform": sys.platform}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="record this environment as the reference")
    a = ap.parse_args(argv)

    F = features()
    digest = hashlib.sha256(np.ascontiguousarray(F, dtype=np.float32)
                            .tobytes()).hexdigest()
    v = versions()

    if a.write:
        json.dump({
            "what": "SHA256 of the float32 (%d, 1038) ligand feature matrix the "
                    "model is fed, for the %d SMILES in familyfm/rdkit_check.py"
                    % (len(SMILES), len(SMILES)),
            "encoder": "Morgan COUNT fingerprint radius 2, 1024 bits via "
                       "rdFingerprintGenerator.GetMorganGenerator, then 14 "
                       "descriptors: MolWt, HeavyAtomCount, NumBonds, "
                       "NumRotatableBonds, RingCount, then counts of "
                       "C N O S F Cl Br I P",
            "shape": list(F.shape),
            "sha256": digest,
            "row_sums": [round(float(x), 4) for x in F.sum(axis=1)],
            "nonzero_per_row": [int(x) for x in (F != 0).sum(axis=1)],
            "reference_environment": v,
            "how_to_verify": "python -m familyfm.rdkit_check",
        }, open(DEST, "w"), indent=1)
        print("wrote %s" % DEST)
        print("sha256 %s" % digest)
        for k, val in v.items():
            print("  %-14s %s" % (k, val))
        return 0

    ref = json.load(open(DEST))
    print("reference built under:")
    for k, val in ref["reference_environment"].items():
        print("  %-14s %s" % (k, val))
    print("this environment:")
    for k, val in v.items():
        print("  %-14s %s" % (k, val))
    if digest == ref["sha256"]:
        print("\nMATCH. This environment computes byte-identical ligand features.")
        return 0
    print("\nMISMATCH. sha256 %s, expected %s" % (digest, ref["sha256"]))
    got = [round(float(x), 4) for x in F.sum(axis=1)]
    bad = [i for i, (g, r) in enumerate(zip(got, ref["row_sums"])) if g != r]
    if bad:
        print("rows that differ (index, this, reference):")
        for i in bad:
            print("   %2d  %12.4f  %12.4f  %s" % (i, got[i], ref["row_sums"][i],
                                                  SMILES[i]))
    else:
        print("every row sum agrees, so the difference is below the row-sum "
              "rounding: a float or a single-bit difference, not a version gap.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
