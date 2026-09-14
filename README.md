# Family Foundation Model

**Two pairwise comparators over a roster spanning 38 protein families across
the two models.**
Neither predicts a value. Each answers one two-class question and returns a
preference probability.

![The two models side by side: preference on the left, one target and two
compounds; target preference on the right, one compound and two targets from
different families.](docs/familyfm-pair.png)

Both answer a question about **preference**, and they differ in whose preference
is being asked about.

| model | you supply | it answers | row layout | width |
|---|---|---|---|---|
| **compound preference, "LSL"** | one target, two compounds | which compound that target prefers | ligand, sequence, ligand | 2,556 |
| **target preference, "SLS"** | one compound, two targets from different families | which target that compound prefers | sequence, ligand, sequence | 1,998 |

Target preference asks where a compound goes; compound preference asks which
compound to take there. Comparing two targets against one compound is what
medicinal chemistry calls **selectivity**, and the word is kept for that
literature rather than used as the name of either model here.

Each answer is one probability, and **the probability is its own confidence**:
the further it sits from 0.5, the more often it turns out right. Every figure
below is printed with the measured accuracy of its strength band.

**Project home: <https://familyfoundationmodel.com>** - both models in a browser
with no install, every target they cover, and the methods behind every number
here.

**Model release: 13 September 2026.** Trained on **ChEMBL 37 alone**, which is what
makes the weights freely downloadable.

## What "prefers" means, and why that word

The endpoints are Ki, Kd, IC50, EC50 and Kb, and they do not all measure one
physical quantity, so "binds more tightly" would be literally true only for Kd
and Ki and "predicts potency" would claim a number this model never returns.

Two properties of the construction make **prefers** exact rather than a hedge.
Every comparison is **endpoint matched**, so the two readings always share one
endpoint and a Ki is never ranked against an IC50. And every reading is converted
to **pActivity**, minus the log of a concentration, so a larger value means a
lower concentration for all five endpoints, whether the reading is an inhibition
endpoint or an agonist one. There is no direction to flip between them.

The operational form: **both readings in a comparison come from the same target
and the same endpoint, and the model predicts which one is active at the lower
concentration.** It ranks a pair. Two compounds can both be nanomolar and the
answer still come back at 0.97.

## Compound preference, "LSL"

Given **one target and two compounds**, which compound that target prefers. This
is the question you ask of a hit list or a congeneric series: rank a purchasable
set before committing anything to assay.

![Two compounds, molibresib and inobrodib, put to one target, the histone
acetyltransferase p300. Each compound becomes a 1,024-bit Morgan count
fingerprint plus 14 descriptors and the target sequence becomes 480 numbers
through ESM2. All three blocks enter a random forest as one row, which returns
which of the two compounds the target prefers.](docs/arch-preference.svg)

**The worked case, and it is a real held-out one.** Histone acetyltransferase p300,
`Q09472`, against Molibresib and Inobrodib. Both are
clinical-stage bromodomain compounds, so this is a hard discrimination inside one
mechanistic class. p300 carries two family labels, reader and writer: its
acetyltransferase domain writes acetyl marks onto lysines and its bromodomain
reads marks already there, so one protein does both jobs. Both readings are
Kd, so the comparison is endpoint matched: pKd
**4.149** against **8.77**. The model answers
**0.790** for Inobrodib, the measured order, and neither
compound appeared anywhere in training.

```
Compound A CCNC(=O)C[C@@H]1N=C(c2ccc(Cl)cc2)c2cc(OC)ccc2-n2c(C)nnc21   Molibresib
Compound B CO[C@H]1CC[C@H](n2c([C@@H]3CCCC(=O)N3c3ccc(F)c(F)c3)nc3cc(-c4c(C)noc4C)ccc32)CC1   Inobrodib
Target     Q09472   EP300, Histone acetyltransferase p300
```

### Multi-family, and what that claims

A single preference comparison never crosses a family boundary, because the
target is held fixed. The model is multi-family because its **roster** is: one
forest ranks compounds at any of 2,079 targets across
36 protein families. That is the larger of the two claims
in this repository.

### The numbers

Holdout is compound-disjoint on **both** compounds by InChIKey, with straddling
comparisons discarded rather than assigned. Every figure is measured on
65,725 held-out comparisons in which neither compound was seen
during fitting.

| | |
|---|---|
| **Accuracy** | **0.710** on 65,725 held-out comparisons over 43,961 compounds |
| Targets servable | 2,079 across 36 protein families |
| Training comparisons | 1,827,578 over 474,708 compounds and 2,077 targets |
| Row layout | ligand, sequence, ligand, 2,556 columns |

### Read the headline beside its baselines

A preference comparator can score respectably while ignoring the target
entirely, because some compounds simply read lower everywhere.

| what is being asked | accuracy |
|---|---|
| **the model** | **0.710** |
| the same forest with the **target removed** | 0.665 |
| always pick the heavier compound, the size trap | 0.569 |
| majority, which the swap fixes at a coin flip | 0.500 |

The obvious compound-history baseline, each compound's mean pActivity across
every target it was measured at, is a **lookup rather than a baseline**: the
average includes the very measurement being predicted, and most held-out
compounds are measured at exactly one target. Taken honestly, from targets
**other than** the one being compared at, it can answer only 25,374 of the
65,725 held-out comparisons, 38.6% of them, because it
needs both compounds measured somewhere else already. The model needs only the
two structures and the target sequence.

| strength gate | comparisons both can answer | the model | the history baseline |
|---|---|---|---|
| answer everything | 25,374 | **0.696** | 0.672 |
| 0.60 | 9,961 | **0.834** | 0.724 |
| 0.70 | 3,771 | **0.914** | 0.784 |
| 0.80 | 1,195 | **0.963** | 0.829 |
| 0.90 | 270 | **0.993** | 0.874 |

### Accuracy rises with prediction strength

![Accuracy rises from 0.710 answering everything to
0.985 at a strength cutoff of 0.90, while the share of
comparisons still answered falls from 100 percent to
1.2%.](docs/strength-tradeoff-preference.svg)

| strength at or above | comparisons kept | share of the held-out set | accuracy |
|---|---|---|---|
| answer everything | 65,725 | 100.0% | 0.710 |
| 0.60 | 28,097 | 42.8% | 0.844 |
| 0.70 | 11,435 | 17.4% | 0.919 |
| 0.80 | 3,995 | 6.1% | 0.961 |
| 0.90 | 794 | 1.2% | 0.985 |

**0.80 is the sensible default.** Show strength beside every prediction and read
it with the accuracy of its band, never on its own.

### Every family with held-out comparisons

The roster this model can score spans 36 families; 31 of them carry enough
held-out comparisons to quote an accuracy for, and those 31 are the table below,
summing to all 65,725 held-out comparisons. This is what evidences the
multi-family claim. The holdout fraction is set per
family, aiming each family at a comparable absolute number of held-out
comparisons rather than the same share of its own pool: families differ in pool
size by a factor of 600, so equal shares would have left the smallest with too
few to quote. Every accuracy is printed with its n. The same table at every gate
is on the site, with a control that moves the gate:
<https://familyfoundationmodel.com/preference.html>

| family | held out | accuracy, no gate | n at 0.80 | accuracy at 0.80 |
|---|---|---|---|---|
| Kinase | 9,585 | 0.688 | 497 | 0.958 |
| Family A G protein-coupled receptor | 7,266 | 0.719 | 485 | 0.930 |
| Protease | 3,160 | 0.712 | 276 | 0.953 |
| Transferase | 2,517 | 0.712 | 187 | 0.952 |
| Electrochemical transporter | 2,306 | 0.668 | 105 | 0.924 |
| Hydrolase | 2,069 | 0.700 | 116 | 0.948 |
| Fatty acid binding protein family | 2,019 | 0.687 | - | - |
| Reader | 2,012 | 0.719 | 190 | 0.958 |
| Other ion channel | 1,993 | 0.752 | 243 | 0.979 |
| Ligand-gated ion channel | 1,977 | 0.686 | 77 | 0.974 |
| Eraser | 1,974 | 0.720 | 88 | 0.966 |
| Ligase | 1,952 | 0.732 | 95 | 0.968 |
| Phosphatase | 1,944 | 0.749 | 118 | 1.000 |
| Writer | 1,931 | 0.795 | 174 | 0.989 |
| Oxidoreductase | 1,918 | 0.721 | 157 | 0.955 |
| Isomerase | 1,914 | 0.732 | 217 | 0.977 |
| Lyase | 1,893 | 0.730 | 110 | 0.991 |
| Toll-like and Il-1 receptors | 1,880 | 0.695 | 30 | 1.000 |
| Nuclear receptor | 1,871 | 0.718 | 89 | 0.989 |
| Phosphodiesterase | 1,829 | 0.718 | 76 | 0.987 |
| Cytochrome P450 | 1,785 | 0.696 | - | - |
| Primary active transporter | 1,781 | 0.685 | 75 | 1.000 |
| Voltage-gated ion channel | 1,773 | 0.714 | 110 | 0.973 |
| Family C G protein-coupled receptor | 1,751 | 0.697 | 83 | 0.988 |
| Family B G protein-coupled receptor | 1,727 | 0.745 | 316 | 0.937 |
| Aminoacyltransferase | 960 | 0.663 | - | - |
| Frizzled family G protein-coupled receptor | 515 | 0.742 | - | - |
| Taste family G protein-coupled receptor | 474 | 0.593 | - | - |
| Calcium channel auxiliary subunit alpha2delta family | 434 | 0.624 | - | - |
| Transmembrane 1-electron transfer carriers | 367 | 0.730 | - | - |
| Group translocator | 148 | 0.723 | - | - |

### By endpoint

The preference model carries real EC50 and Kb volume, so it reaches agonist and
functional readings as well as binding ones.

| endpoint | comparisons | accuracy |
|---|---|---|
| IC50 | 34,481 | 0.713 |
| Ki | 12,900 | 0.717 |
| EC50 | 9,375 | 0.713 |
| Kd | 8,431 | 0.685 |
| Kb | 538 | 0.731 |

## Target preference, "SLS"

Given **one ligand and two protein targets from different families**, which
target the ligand prefers.

![Olanzapine put to two targets from different protein families: the histamine
H1 receptor, drawn as a seven-helix bundle in a membrane, and hERG, drawn as four
subunits around a central pore. Each target sequence becomes 480 numbers through
ESM2 and the ligand becomes a 1,024-bit Morgan count fingerprint plus 14
descriptors. All three blocks enter a random forest as one row, which returns
which of the two targets the compound prefers.](docs/arch-cross-family.svg)

**The worked case.** Olanzapine put to the histamine H1 receptor, a G
protein-coupled receptor, and to hERG, a voltage-gated ion channel. Both are real
measurements from the training data: pKi **8.50** at H1 against **4.44** at hERG,
so the compound is active at a concentration 4.1 log units lower at H1. The model answers **0.902 for
H1**, the measured order. The target pictures are schematics of each protein's
architecture, not structures.

```
SMILES   Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1      olanzapine
Target A P35367   HRH1,  Family A G protein-coupled receptor
Target B Q12809   KCNH2, Voltage-gated ion channel
```

Paste that into <https://familyfoundationmodel.com/rank.html> and you should get
the same answer.

The whole comparison is one row and **the order is the question**. Nothing is
scored on its own and then subtracted.

### What target preference is for

**Off-target triage and repurposing.** Give it a compound and a panel of targets
drawn from different families, and it ranks which the compound leans toward.

It is **not** a safety or toxicity screen. Breadth across families in this data
partly records how many assay panels a compound went through rather than how
promiscuous it is, and reading it as toxicity would be a claim the measurements
do not support.

### The target-preference numbers

Compound-disjoint holdout by InChIKey: **zero ligands appear on both sides**, so
every held-out comparison involves chemistry the model was never fitted on.

| | |
|---|---|
| **Accuracy** | **0.750** on 8,689 held-out comparisons over 2,195 ligands |
| Targets servable | 1,879 across 38 protein families |
| Comparisons, total | 89,888 over 22,588 ligands and 290 family pairings |
| Comparisons fitted on | 81,199 over 20,393 ligands |

#### Read the headline beside its baselines. This is the point.

A cross-family comparator can look respectable while knowing nothing about the
compound, because some families simply carry more potent chemistry than others.

| what is being asked | accuracy |
|---|---|
| **the model** | **0.750** |
| always pick whichever family usually wins that pairing | 0.655 |
| the same forest with the **ligand removed** | 0.710 |

The compound contributes about **four points** on top of target identity. That
margin was reproduced at 3.9, 4.5 and 4.1 points across three independently built
versions of this model. It is real, and it is smaller than the headline suggests
on its own.

#### Accuracy rises with prediction strength

![Accuracy rises from 0.750 answering everything to 0.966 at a strength cutoff of
0.90, while the share of comparisons still answered falls from 100 percent to 22
percent.](docs/strength-tradeoff.svg)

Strength is the larger of the two returned probabilities, so it runs 0.5 to 1.0.

| strength at or above | comparisons kept | accuracy |
|---|---|---|
| answer everything | 100% | 0.750 |
| 0.60 | 74.4% | 0.815 |
| 0.70 | 52.0% | 0.877 |
| 0.80 | 35.9% | 0.933 |
| 0.90 | 22.1% | 0.966 |

#### A near-tie is a near-tie

| true separation | comparisons | accuracy |
|---|---|---|
| under half a log | 2,283 | 0.573 |
| half a log to one log | 1,766 | 0.701 |
| one to two logs | 2,359 | 0.797 |
| beyond two logs | 2,281 | 0.915 |

#### Some pairings are much harder than others

Across pairings carrying at least 100 held-out comparisons, strongest first:

| pairing | comparisons | accuracy |
|---|---|---|
| Kinase against Voltage-gated ion channel | 144 | 0.931 |
| Cytochrome P450 against Kinase | 278 | 0.899 |
| Cytochrome P450 against Oxidoreductase | 105 | 0.886 |

and the hardest:

| pairing | comparisons | accuracy |
|---|---|---|
| Electrochemical transporter against Family A G protein-coupled receptor | 531 | 0.655 |
| Cytochrome P450 against Voltage-gated ion channel | 121 | 0.744 |
| Eraser against Kinase | 121 | 0.744 |

Pairings with fewer than 30 held-out comparisons are not quoted, here or in the
manifest.

#### By endpoint, target preference

| endpoint | comparisons | accuracy |
|---|---|---|
| IC50 | 6,577 | 0.776 |
| Ki | 1,440 | 0.697 |
| Kd | 558 | 0.588 |
| EC50 | 114 | 0.693 |

The Kd row is thin and weak; do not lean on it.

## Install

```bash
git clone https://github.com/smuskal/ffm.git
cd ffm
./install.sh                          # both bundles
FFM_MODELS=selectivity ./install.sh   # or just the smaller one, the target-preference bundle
FFM_MODELS=preference  ./install.sh
```

**Nothing is fetched from Hugging Face, and nothing needs torch.** Each bundle
ships `sequence_vectors.npz`, the ESM2 vectors for every target it can score,
already computed. `predict.py` reads them out of that file; asked about an
accession the bundle does not carry it raises rather than reaching for a
checkpoint, so scoring is fully offline and cannot be broken by an upstream
model being moved or relicensed. `embed.py`, `torch` and `transformers` are
needed only to BUILD a bundle for a sequence that is not already in one, which
is why they are commented out of `requirements.txt`.

That builds the environment with the versions pinned below, fetches the bundles,
and **proves each one works** by checking its forest against the checksum in its
own manifest and replaying the reference predictions that ship inside it. If they
do not reproduce to 1e-6 the install fails rather than reporting success.

Hardware: the target-preference forest is **0.09 GB on disk** and about **1.0 GB
resident**, loading in under four seconds, so an 8 GB laptop is comfortable. The
preference forest is **1.18 GB on disk** and loads to rather more than that, so
give it 16 GB. Load either **once** at process start, never per request. No GPU
is used.

## Run it

Rank targets for one compound, the target-preference question:

```bash
./ffm.sh rank --smiles "CC(=O)Oc1ccccc1C(=O)O" --targets P00533 P08684 P28223
./ffm.sh compare --smiles "CC(=O)Oc1ccccc1C(=O)O" --a P00533 --b P08684
```

Rank compounds at one target, the compound-preference question:

```bash
./ffm.sh rank-ligands --target P00533 --smiles-file hits.smi
./ffm.sh compare-ligands --target P00533 --a "<SMILES>" --b "<SMILES>"
```

`./ffm.sh` activates the virtual environment `install.sh` built and then runs
`python -m familyfm.cli`, so it works without installing anything system wide.
Equivalently, once that environment is active,
`python -m familyfm.cli rank --smiles ...`.

`hits.smi` is one SMILES per line, with an optional name in a second column.
Each subcommand picks the bundle it needs; `--bundle` overrides that.

**Targets are named by UniProt accession, not gene symbol.** Gene symbols are
many-to-many against accessions and the failure is silent: in an earlier roster
`CHRM5` resolved to an accession belonging to a different protein entirely. The
full roster is at <https://familyfoundationmodel.com/targets.html>.

### Or from Python

```python
import familyfm.bundle_predict as P
m = P.load("path/to/xfam_v1")

P.compare_targets(m, smiles, "P00533", "P08684")   # P(the ligand prefers the first)
P.rank_targets(m, smiles, ["P00533", "P08684", "P28223"])
P.family_of(m, "P00533")                            # 'Kinase'
P.targets(m)                                        # every servable accession
```

The compound-preference bundle loads the same way and exposes the mirror calls:

```python
import familyfm.preference_predict as P
m = P.load("path/to/lsl_v2_stratified")

P.compare_ligands(m, smiles_a, smiles_b, "P00533")   # P(the target prefers the first)
P.compare_many_ligands(m, candidates, reference, "P00533")
P.rank_ligands(m, [s1, s2, s3, s4], "P00533")        # preference order at one target
P.strength(p)                                        # max(p, 1-p)
```

Unparseable SMILES **raise** rather than being silently dropped, so catch the
`ValueError` and tell the user which input failed.

`rank_targets` scores every pair among the targets given, both orders averaged,
and returns them ranked. With n targets that is n(n-1)/2 comparisons, so keep
panels to about a dozen per call.

**Always batch.** `compare_many_targets` amortizes the forest call across
ligands. Load the bundle **once**, never per request.

## Two properties you can rely on

Both are asserted by the bundle's own self-test:

* **Exact antisymmetry.** Every prediction averages both input orders, so
  `compare_targets(A,B) + compare_targets(B,A) == 1` to machine precision, and a
  target compared with itself returns exactly 0.5. Bypassing that average is the
  one bug that silently halves this model's value.
* **Endpoint matching at training time.** An IC50 is only ever compared with
  another IC50, so there is no endpoint argument at inference and none is needed.

## How a target-preference comparison is formed

Four rules, and they are the substance of the method:

1. **Endpoint matched**, so the comparison is between two readings of the same kind.
2. **The two targets share no family label**, which is what makes this cross-family.
3. **Ties are dropped.** Two equal readings order nothing, and entering them
   teaches the forest that whichever side is written first wins.
4. **Which target is A is randomized**, and every comparison is entered twice
   with the sequence blocks exchanged and the answer inverted.

Row layout is sequence, ligand, sequence: 480 + 1,038 + 480 = **1,998** columns.
The ligand block is a Morgan count fingerprint, radius 2, 1,024 bits, plus 14
descriptors. Sequences are ESM2 `esm2_t12_35M_UR50D`, mean pooled over residues.

## Pin these versions

A joblib forest is a pickled object graph: a different scikit-learn either
refuses to load it or loads it and scores differently. RDKit computes 1,024 of
every ligand block's 1,038 dimensions, so a different RDKit changes the
fingerprint and moves every score **with nothing to detect it**.

| package | pinned |
|---|---|
| numpy | 2.2.6 |
| python | 3.10.15 |
| sklearn | 1.7.2 |

`torch` and `transformers` are needed only to embed a sequence the bundle does
not already carry. For the shipped roster the vectors are precomputed, so normal
use never imports them.

## Reproducing the training set

Every step runs from ChEMBL 37 and nothing else, which is the point of building
it this way: anyone with the same release can rebuild it.

```bash
# 1. scope ChEMBL to human single-protein targets carrying a family
python -m familyfm.extract --db data/chembl_37.db --out data/measurements.csv

# 2. form endpoint-matched cross-family comparisons
python -m familyfm.pairs --measurements data/measurements.csv --out data/pairs.csv

# 3. fit, hold out compound-disjoint, measure
python -m familyfm.train --pairs data/pairs.csv --out bundles/xfam

# 4. write the bundle and its reference predictions
python -m familyfm.finalize --bundle bundles/xfam --measurements data/measurements.csv
```

`familyfm/features.py` is imported by the steps above rather than run on its
own: it is the feature recipe, not a stage. Training needs `pandas` as well as
the pinned four; see the training section of `requirements.txt`.

`familyfm/rdkit_check.py` verifies that your RDKit produces the same ligand
features the forest was fitted on, which is the check that catches a silent
scoring regression.

## What is in the bundle

| file | what it is |
|---|---|
| `familyfm_selectivity.joblib` or `familyfm_preference.joblib` | the forest, one per bundle |
| `predict.py` | the inference API; load it by file path, nothing to install |
| `sequence_vectors.npz` | ESM2 vectors for every servable target, precomputed |
| `sequence_index.json` | sequence key to row |
| `target_index.json` | accession to sequence key, and the family per target |
| `MANIFEST.json` | provenance, hyperparameters, measured performance |
| `reference_predictions.json` | 30 fixed inputs and the outputs this bundle produced |
| `heldout_scored.csv` | every held-out comparison with its probability |

Because `heldout_scored.csv` ships with the bundle, **every number above can be
recomputed from the download** without rerunning the model.

## Data

ChEMBL 37 only: 1,263,626 measurements over
2,705 targets and 794,638 ligands
after scoping to human single-protein targets carrying a protein family.
Endpoints are Ki, IC50, Kd, EC50 and Kb.

Scope is curated dose-response measurement: Ki, IC50, Kd, EC50 and Kb read from
assays reporting a concentration, taken from the medicinal-chemistry literature
and patents that ChEMBL curates. The effect of holding to that scope was
measured, and the numbers behind it are published at
<https://familyfoundationmodel.com/methods.html>.

Activity is derived from the reported value and unit rather than read from a
precomputed column. Duplicates collapse by **median**, never by the most potent
value, because the best value selects unit errors. Censored readings are dropped.

## Licensing

| | |
|---|---|
| code in this repository | **Apache 2.0** (`LICENSE`) |
| model weights, both bundles | **CC BY-SA 4.0** (`LICENSE-MODEL.txt`) |

The weights carry share-alike because ChEMBL does. This model is trained on
ChEMBL 37 and nothing else, so the weights are a derived work of a CC BY-SA
dataset and pass those conditions forward. Using the **code** with your own data
puts no such obligation on you.

The sibling [Kinase Foundation Model](https://kinasefoundationmodel.com/) ships
its weights under research-and-evaluation-only terms because it is trained on a
commercial database. This one is not, which is why it can be downloaded openly.

## Related work

* [Kinase Foundation Model](https://kinasefoundationmodel.com/) - the same two
  layouts within protein kinases
* [GPCR Foundation Model](https://gpcrfoundationmodel.com/) - the same two
  layouts within G protein-coupled receptors
* [PharmCast](https://pharmcast.ai/) and
  [Reverse Screen](https://reversescreen.ai/)

## Citing

> Family Foundation Model, Eidogen-Sertanty, Inc., 2026. Model release
> 13 September 2026. Trained on ChEMBL 37.
> https://familyfoundationmodel.com

Cite ChEMBL alongside it.

(c) 2026 Eidogen-Sertanty, Inc.
