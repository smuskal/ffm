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
#   FFM_MODELS=selectivity ./install.sh     the target comparator, 0.42 GB
#   FFM_MODELS=preference  ./install.sh     the compound comparator, 1.18 GB
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
# The bundles are served from the project site. Nothing here names any hosting
# provider or account: the site is the published address and stays the published
# address even if what sits behind it changes.
# Downloads go through the site's download endpoint, which records them and
# redirects to the archive, so -L is required. Set FFM_BASE_URL to the static
# /models path to bypass it.
BASE_URL="${FFM_BASE_URL:-https://familyfoundationmodel.com/api/v1/download}"
WANT="${FFM_MODELS:-both}"

# Python 3.12 is the recommended interpreter; 3.10 to 3.13 work. The pins in
# requirements.txt have wheels for exactly that range on macOS and Linux.
# Set FFM_PYTHON to choose one yourself.
pick_python () {
  local c
  for c in ${FFM_PYTHON:-} python3.12 python3.13 python3.11 python3.10 python3; do
    command -v "$c" >/dev/null 2>&1 || continue
    if "$c" -c 'import platform, sys
import ensurepip, venv
sys.exit(0 if (3, 10) <= sys.version_info[:2] <= (3, 13) else 1)' >/dev/null 2>&1
    then command -v "$c"; return 0; fi
    if [ -n "${FFM_PYTHON:-}" ]; then return 1; fi
  done
  return 1
}
if ! PY="$(pick_python)"; then
  cat >&2 <<'MSG'
No suitable Python found. Install needs Python 3.10 to 3.13; 3.12 is a good choice.
  macOS:  brew install python@3.12      then rerun ./install.sh
  Linux:  your package manager's python3.12, or  uv python install 3.12
  any:    FFM_PYTHON=/path/to/python3.12 ./install.sh
MSG
  exit 1
fi

echo "== environment"
echo "   using $PY ($("$PY" -c 'import platform,sys;print(sys.version.split()[0], platform.machine())'))"
# a .venv left by an earlier run under a different Python is rebuilt, not reused
if [ -x "$HERE/.venv/bin/python" ] && \
   [ "$("$HERE/.venv/bin/python" -c 'import sys;print(sys.version_info[:2])')" != "$("$PY" -c 'import sys;print(sys.version_info[:2])')" ]; then
  rm -rf "$HERE/.venv"
fi
"$PY" -m venv "$HERE/.venv"
# shellcheck disable=SC1091
. "$HERE/.venv/bin/activate"
pip install --quiet --upgrade "pip>=24.1"
if ! pip install --quiet -r "$HERE/requirements.txt"; then
  echo >&2
  echo "pip could not install the pinned libraries. What this environment is:" >&2
  python -VV >&2
  python -c 'import platform, sys; print("executable", sys.executable); print("machine   ", platform.machine())' >&2
  pip --version >&2
  echo "Upgrade pip first (pip install -U pip) and rerun. On a Mac, machine x86_64" >&2
  echo "means an Intel Python; a native arm64 Python is the usual fix." >&2
  exit 1
fi

# The reference replay below cannot see every encoder change, so the ligand
# features are checked directly, before anything is downloaded.
echo "== proving RDKit computes the features the models were built on"
if ! python -m familyfm.rdkit_check > "$HERE/.rdkit_check.log" 2>&1; then
  cat "$HERE/.rdkit_check.log" >&2
  echo "The ligand features differ from the build. Reinstall with the pins in" >&2
  echo "requirements.txt rather than an RDKit already on this machine." >&2
  exit 1
fi
echo "   features    OK  byte-identical, rdkit $(python -c 'import rdkit; print(rdkit.__version__)')"

fetch () {          # fetch <key> <archive-file> <bundle-dir> <joblib>
  local key="$1" file="$2" dir="$HERE/ffm-models/$3" jb="$4"
  # the download endpoint takes the key; the static /models path takes the file name
  local url="$BASE_URL/$key"
  case "$BASE_URL" in */models|*/models/) url="${BASE_URL%/}/$file" ;; esac
  if [ -f "$dir/$jb" ]; then
    echo "   $3 already present, not refetched"
    return 0
  fi
  mkdir -p "$HERE/ffm-models"
  local tmp; tmp="$(mktemp "${TMPDIR:-/tmp}/ffm.XXXXXX")"
  # -L follows the redirect the site issues; --retry covers a dropped connection
  # part way through what can be a gigabyte
  curl -fL --retry 3 "$url" -o "$tmp"
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
  fetch target xfam_v2.tar.gz xfam_v2_samefamily familyfm_selectivity.joblib
  echo "== proving the target-preference model works"
  verify xfam_v2_samefamily familyfm_selectivity.joblib selectivity
fi

if [ "$WANT" = "both" ] || [ "$WANT" = "preference" ]; then
  echo "== compound-preference bundle"
  fetch preference lsl_v1.tar.gz lsl_v2_stratified familyfm_preference.joblib
  echo "== proving the compound-preference model works"
  verify lsl_v2_stratified familyfm_preference.joblib preference
fi

echo
echo "Installed. Try:"
echo "  . .venv/bin/activate"
echo "  ./ffm.sh rank --smiles 'CC(=O)Oc1ccccc1C(=O)O' --targets P00533 P08684 P28223"
echo "  ./ffm.sh rank-ligands --target P00533 --smiles-file hits.smi"
