"""The target-preference inference API, re-exported under its historical name.

`predict.py` is the file that ships INSIDE a bundle, where it is loaded by path
with nothing installed. This module is the same API under the name the README and
the CLI have always used. It is a re-export rather than a copy, because the two
were byte-identical copies and a copy drifts silently.
"""
from familyfm.predict import (           # noqa: F401
    load, targets, family_of, ligand_features, seq_vector, strength,
    compare_targets, compare_many_targets, rank_targets,
)
