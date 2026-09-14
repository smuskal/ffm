#!/usr/bin/env bash
# Install the Family Foundation Model: environment, code and both model bundles.
#
# It PROVES the install rather than reporting success: each bundle ships its own
# reference predictions and this replays them. If they do not reproduce to the
# bundle's tolerance the install fails, because a forest loaded under the wrong
# scikit-learn or scored with the wrong RDKit still answers, just differently,
# and nothing else would notice.
#
# Set FFM_MODELS to install only one of them:
#   FFM_MODELS=selectivity ./install.sh     the smaller forest, 0.09 GB
#   FFM_MODELS=preference  ./install.sh     the larger forest, 1.18 GB
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
# The bundles are served from the project site. Nothing here names any hosting
# provider or account: the site is the published address and stays the published
# address even if what sits behind it changes.
BASE_URL="${FFM_BASE_URL:-https://familyfoundationmodel.com/models}"
WANT="${FFM_MODELS:-both}"

echo "== environment"
python3 -m venv "$HERE/.venv"
# shellcheck disable=SC1091
. "$HERE/.venv/bin/activate"
pip install --quiet --upgrade "pip>=24.1"
pip install --quiet -r "$HERE/requirements.txt"

fetch () {          # fetch <archive> <bundle-dir> <joblib>
  local archive="$1" dir="$HERE/ffm-models/$2" jb="$3"
  if [ -f "$dir/$jb" ]; then
    echo "   $2 already present, not refetched"
    return 0
  fi
  mkdir -p "$HERE/ffm-models"
  local tmp; tmp="$(mktemp -t ffm)"
  # -L follows the redirect the site issues; --retry covers a dropped connection
  # part way through what can be a gigabyte
  curl -fL --retry 3 "$BASE_URL/$archive" -o "$tmp"
  tar xzf "$tmp" -C "$HERE/ffm-models"
  rm -f "$tmp"
}

verify () {         # verify <bundle-dir> <joblib> <selectivity|preference>
  python3 - "$HERE/ffm-models/$1" "$2" "$3" <<'PY'
import hashlib, importlib.util, json, sys
B, JOBLIB, KIND = sys.argv[1], sys.argv[2], sys.argv[3]

man = json.load(open(B + "/MANIFEST.json"))
h = hashlib.sha256()
with open(B + "/" + JOBLIB, "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 20), b""):
        h.update(chunk)
if h.hexdigest() != man["model_sha256"]:
    sys.exit("checksum mismatch: the download is not the model this manifest describes")
print("   sha256      OK")

spec = importlib.util.spec_from_file_location("fp_" + KIND, B + "/predict.py")
fp = importlib.util.module_from_spec(spec); sys.modules["fp_" + KIND] = fp
sys.path.insert(0, B); spec.loader.exec_module(fp)
m = fp.load(B)

ref = json.load(open(B + "/reference_predictions.json"))
if KIND == "selectivity":
    call = lambda r: fp.compare_targets(m, r["smiles"], r["target_a"], r["target_b"])
    first = ref["predictions"][0]
    same = fp.compare_targets(m, first["smiles"], first["target_a"], first["target_a"])
    ab = fp.compare_targets(m, first["smiles"], first["target_a"], first["target_b"])
    ba = fp.compare_targets(m, first["smiles"], first["target_b"], first["target_a"])
else:
    call = lambda r: fp.compare_ligands(m, r["smiles_a"], r["smiles_b"], r["target"])
    first = ref["predictions"][0]
    same = fp.compare_ligands(m, first["smiles_a"], first["smiles_a"], first["target"])
    ab = fp.compare_ligands(m, first["smiles_a"], first["smiles_b"], first["target"])
    ba = fp.compare_ligands(m, first["smiles_b"], first["smiles_a"], first["target"])

worst = max(abs(call(r) - r["p"]) for r in ref["predictions"])
if worst > ref["tolerance"]:
    sys.exit("reference predictions differ by %.3e, tolerance %g. Check your "
             "scikit-learn and RDKit versions against requirements.txt."
             % (worst, ref["tolerance"]))
print("   reference   OK  %d predictions, worst %.3e" % (len(ref["predictions"]), worst))

if abs(same - 0.5) > 1e-12 or abs(ab + ba - 1.0) > 1e-12:
    sys.exit("the both-orders average is broken")
print("   symmetry    OK  self %.6f, A+B %.6f" % (same, ab + ba))
print("   targets     %d servable" % len(fp.targets(m)))
PY
}

if [ "$WANT" = "both" ] || [ "$WANT" = "selectivity" ]; then
  echo "== target-preference bundle"
  fetch xfam_v1.tar.gz xfam_v1 familyfm_selectivity.joblib
  echo "== proving the target-preference model works"
  verify xfam_v1 familyfm_selectivity.joblib selectivity
fi

if [ "$WANT" = "both" ] || [ "$WANT" = "preference" ]; then
  echo "== compound-preference bundle"
  fetch lsl_v1.tar.gz lsl_v2_stratified familyfm_preference.joblib
  echo "== proving the compound-preference model works"
  verify lsl_v2_stratified familyfm_preference.joblib preference
fi

echo
echo "Installed. Try:"
echo "  . .venv/bin/activate"
echo "  ./ffm.sh rank --smiles 'CC(=O)Oc1ccccc1C(=O)O' --targets P00533 P08684 P28223"
echo "  ./ffm.sh rank-ligands --target P00533 --smiles-file hits.smi"
