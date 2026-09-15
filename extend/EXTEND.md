# Adding your own data to a released Family Foundation Model

`ffm_extend.py` fits trees on measurements you hold and merges them into a
released model. It runs entirely on your machine. Your structures, sequences and
activity values are never written into the bundle and never leave your computer;
only trees fitted from them are.

Both released models are supported, and the tool reads the layout from the
bundle's own manifest:

| | compound preference | target preference |
|---|---|---|
| bundle | `lsl_v2_stratified` | `xfam_v1` |
| layout | LSL, `[ ligand A \| sequence \| ligand B ]` | SLS, `[ sequence A \| ligand \| sequence B ]` |
| the question | which of two compounds a target prefers | which of two targets a compound prefers |

    python ffm_extend.py --base ./lsl_v2_stratified --data mine.csv --out ./mine_lsl
    python ffm_extend.py --base ./xfam_v1           --data mine.csv --out ./mine_sls

The base bundle is opened read-only and is never modified. The output is a
complete bundle: it carries the merged forest under the base model's own
filename, so the bundle's `predict.py` loads it with no code change.

---

# Your targets do not have to be on the roster

The released models serve 2,079 and 1,879 targets. You are not limited to those.

**Name a target by `accession`** to use one the bundle already carries.
**Supply its `sequence`** to add a target the released model has never seen. The
sequence is embedded with the same ESM2 recipe the models were built on, your
trees are fitted on it, and the new target is written into the output bundle's
`sequence_vectors.npz` and `target_index.json`. After that the bundle's own
`predict.py` can be asked about it like any other target.

Give the new target an `accession` of your choosing alongside the sequence and it
will be known by that name. Add a `family` column to label it. Neither is
required; without an accession the tool assigns one.

Adding a target needs `torch` and `transformers` installed, because the sequence
has to be embedded. Naming targets only by accession does not.

    pip install torch transformers

---

# Input format

CSV, UTF-8, header required. Column order does not matter and unknown columns are
ignored. Two shapes are accepted and the tool detects which you have.

## Shape 1: one measurement per row

The natural export from a knowledgebase or an assay database, and usually what
you already have.

    smiles,accession,pic50,relation,assay_measure

The tool forms the pairs for you. For compound preference it pairs the compounds
measured on each target; for target preference it pairs the targets each compound
was measured against. The same file builds either model.

A measurement is only usable if something else shares its group: another compound
on the same target for compound preference, another target for the same compound
for target preference. The run reports how many fell out, because a large file can
be almost entirely unpairable for one of the two layouts.

## Shape 2: comparisons you have already paired

**Compound preference (LSL):** two compounds, one target.

    smiles_a,smiles_b,accession,pic50_a,pic50_b,relation_a,relation_b,assay_measure

**Target preference (SLS):** one compound, two targets.

    smiles,accession_a,accession_b,pic50_a,pic50_b,relation_a,relation_b,assay_measure

The released target preference model was fitted and measured on pairs of targets
from different families. Pairs you add from within one family are accepted and
fitted like any other; measure them on your own held-out comparisons, since the
released figures describe cross-family pairs.

Replace any `accession` with a `sequence` column to add a target off the roster.
`sequence_a` and `sequence_b` work the same way for target preference, and you can
mix them: a new target on one side and a roster accession on the other is the
ordinary case for a program asking whether its compound is selective.

## The activity columns

`pic50_a` and `pic50_b` are pActivity, the negative log of a concentration, so a
larger value means active at a lower concentration. Both readings in a comparison
must come from the same endpoint.

`relation_a` and `relation_b` describe the activity: `=` by default, `>` for at
least this active, `<` for at most. A comparison is used only when the two
intervals are disjoint, so a censored reading that carries no order is dropped
rather than guessed. Two identical exact values are a tie, kept and entered both
ways so the model returns about 0.5.

If you have a winner rather than values, supply a `winner` column of `A` or `B`
instead of the two activity columns.

Every usable comparison is entered twice with the label reversed, which is how
both released models were fitted. Do not supply both orders yourself.

---

# Files in `examples/`

Real ChEMBL data, drawn from each model's own held-out comparisons, except the two
templates.

| file | shape | targets |
|---|---|---|
| `extend_lsl_measurements.csv` | one measurement per row | roster accessions |
| `extend_lsl_comparisons.csv` | paired | roster accessions |
| `extend_lsl_new_target_TEMPLATE.csv` | paired | one target off the roster, by sequence |
| `extend_sls_measurements.csv` | one measurement per row | roster accessions |
| `extend_sls_comparisons.csv` | paired | roster accessions |
| `extend_sls_new_target_TEMPLATE.csv` | paired | a new target against roster accessions |
| `breadth_reference_lsl.csv` | paired | 1,500 held-out comparisons over 967 targets |
| `breadth_reference_sls.csv` | paired | 1,200 held-out comparisons over 1,024 target pairs |

**The two TEMPLATE files carry real structures and a real UniProt sequence, and
placeholder activity values.** The sequence is human ubiquitin carboxyl-terminal
hydrolase 27, which is absent from both rosters, so the files demonstrate adding a
target the models have never seen. Replace the activity values with your own; they
are there to show the columns, not to report a measurement.

The two breadth reference files are held-out comparisons spread across the roster.
They are what `--sweep` uses to measure what your added trees cost on targets you
hold no data for.

---

# Choosing how many trees to add

Your trees vote on every prediction the model makes, including targets you hold no
data for. More of them means a larger gain where your data is and a larger risk
everywhere else. `--sweep` measures both instead of guessing.

    python ffm_extend.py --base ./lsl_v2_stratified \
        --data mine_fit.csv --holdout mine_held.csv \
        --out ./mine_lsl --sweep

It fits once at the largest tree count and evaluates prefixes of the same trees,
which is exact, then recommends the largest count whose loss on the breadth
reference stays within `--max-breadth-loss`, 0.005 by default.

## A measured run

1,858 comparisons across five targets, held back from the released compound
preference model, with 797 more from the same targets as the holdout:

| trees added | your vote share | your holdout | gain | breadth reference | cost |
|---|---|---|---|---|---|
| 5 | 1.6% | 0.689 | +0.004 | 0.691 | -0.001 |
| 10 | 3.2% | 0.693 | +0.008 | 0.691 | -0.001 |
| 20 | 6.2% | 0.695 | +0.010 | 0.692 | 0.000 |
| 40 | 11.8% | 0.709 | +0.024 | 0.691 | -0.001 |
| 80 | 21.1% | 0.727 | +0.041 | 0.692 | 0.000 |

The released model reads 0.685 on that holdout and 0.692 on the breadth
reference. Eighty added trees bought 0.041 on the contributor's own targets with
no measured loss across the roster.

The same procedure on the target preference model, 6,000 contributor comparisons
and a 1,200 comparison holdout, moved 0.753 to 0.763 with 40 trees added.

**Both runs used held-out public ChEMBL data.** They show that the
mechanism works and what the trade-off looks like. What your own data buys depends
on your data.

---

# What extension does to the released figures

**The released accuracy and the released strength bands do not describe an
extended model.** The tool withdraws them: the performance block in the output
manifest is replaced with a withdrawal notice, the released figures are kept
beside it for reference only, and `reference_predictions.json` is renamed to
`.base`, because added trees change every prediction.

Measure an extended model on your own held-out comparisons before quoting any
number from it, and re-measure the prediction strength operating point, since the
accuracy attached to each strength band was measured on the released forest.

The output manifest records the extension: when it ran, the layout, how many trees
were added, how many rows were read and used, what was dropped and why, your vote
share, the base model's checksum, any targets added to the roster, and the
scikit-learn version. A bundle therefore carries its own history.

---

# Two constraints

**scikit-learn must match the bundle.** Trees from different minor versions must
not be merged, and the tool refuses rather than producing a model that loads but
misbehaves. The pinned version is in the bundle manifest.

**A thousand usable comparisons is the floor.** Fewer cannot support a useful tree,
and the result would be noise voting on every prediction. `--allow-small` overrides
this for a format check and says so in the log.
