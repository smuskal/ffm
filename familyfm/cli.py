"""Command line for the Family Foundation Model, both comparators.

Rank targets for one compound, the target-preference question:

    ffm rank --smiles "<SMILES>" --targets P00533 P08684 P28223
    ffm compare --smiles "<SMILES>" --a P00533 --b P08684

Rank compounds at one target, the preference question:

    ffm rank-ligands --target P00533 --smiles "<SMILES>" "<SMILES>" "<SMILES>"
    ffm rank-ligands --target P00533 --smiles-file hits.smi
    ffm compare-ligands --target P00533 --a "<SMILES>" --b "<SMILES>"

    ffm targets [--family Kinase]

Targets are UniProt accessions. Gene symbols are many-to-many against accessions
and the failure is silent, so this tool takes accessions only and prints the
family beside each so you can see what you asked for.
"""
import argparse, importlib.util, json, os, sys

_MODELS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "ffm-models")
DEFAULT_BUNDLE = os.environ.get("FFM_BUNDLE", os.path.join(_MODELS, "xfam_v1"))
# the preference comparator is a separate bundle with its own forest and its own
# measured table; the two are never mixed
DEFAULT_PREFERENCE_BUNDLE = os.environ.get(
    "FFM_PREFERENCE_BUNDLE", os.path.join(_MODELS, "lsl_v2_stratified"))
LIGAND_CMDS = ("rank-ligands", "compare-ligands")


def load_bundle(path):
    if not os.path.isdir(path):
        sys.exit("No model bundle at %s.\nRun ./install.sh, or set FFM_BUNDLE to "
                 "a bundle directory." % path)
    name = "ffm_predict_" + os.path.basename(path.rstrip("/"))
    spec = importlib.util.spec_from_file_location(name,
                                                  os.path.join(path, "predict.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    sys.path.insert(0, path)
    spec.loader.exec_module(mod)
    return mod, mod.load(path)


def band(strength, cumulative):
    """The band a strength falls in, and that band's measured accuracy. A 0.52
    and a 0.95 are not the same answer and must not print the same."""
    best = None
    for cut in sorted(cumulative, reverse=True):
        if strength >= float(cut):
            best = (float(cut), cumulative[cut]["accuracy"])
            break
    return best or (0.5, cumulative["0.50"]["accuracy"])


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ffm", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", default=None,
                    help="bundle directory; defaults to the one the subcommand needs")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("rank", help="rank a panel of targets for one compound")
    r.add_argument("--smiles", required=True)
    r.add_argument("--targets", nargs="+", required=True)

    c = sub.add_parser("compare", help="compare exactly two targets")
    c.add_argument("--smiles", required=True)
    c.add_argument("--a", required=True)
    c.add_argument("--b", required=True)

    rl = sub.add_parser("rank-ligands",
                        help="rank compounds at one target, the preference question")
    rl.add_argument("--target", required=True)
    rl.add_argument("--smiles", nargs="*", default=[])
    rl.add_argument("--smiles-file", default=None,
                    help="one SMILES per line; a second column is taken as a name")

    cl = sub.add_parser("compare-ligands", help="compare exactly two compounds at one target")
    cl.add_argument("--target", required=True)
    cl.add_argument("--a", required=True)
    cl.add_argument("--b", required=True)

    t = sub.add_parser("targets", help="list the targets this model can score")
    t.add_argument("--family", default=None)

    a = ap.parse_args(argv)
    bundle = a.bundle or (DEFAULT_PREFERENCE_BUNDLE if a.cmd in LIGAND_CMDS
                          else DEFAULT_BUNDLE)
    a.bundle = bundle
    P, m = load_bundle(bundle)
    perf = json.load(open(os.path.join(a.bundle, "MANIFEST.json")))["MEASURED_PERFORMANCE"]
    cum = perf["strength_cumulative"]

    if a.cmd in LIGAND_CMDS:
        if a.target not in set(P.targets(m)):
            sys.exit("Not in this model: %s\nTargets are UniProt accessions. "
                     "Run `ffm targets --bundle %s` to list them."
                     % (a.target, a.bundle))
        fam = (P.family_of(m, a.target) or "").replace("|", " and ")

        if a.cmd == "compare-ligands":
            p = P.compare_ligands(m, a.a, a.b, a.target)
            st = max(p, 1 - p)
            _cut, bacc = band(st, cum)
            print("%s (%s)" % (a.target, fam))
            print("  prefers the %s compound at strength %.2f"
                  % ("first" if p >= 0.5 else "second", st))
            print("  that band is right %.2f of the time" % bacc)
            return 0

        smiles, names = list(a.smiles), {}
        if a.smiles_file:
            for line in open(a.smiles_file):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(None, 1)
                smiles.append(parts[0])
                if len(parts) > 1:
                    names[parts[0]] = parts[1]
        if len(smiles) < 2:
            sys.exit("Give at least two compounds, with --smiles or --smiles-file.")
        # every pair is scored, both orders averaged, so the cost is quadratic
        pairs = len(smiles) * (len(smiles) - 1) // 2
        print("%d compounds, %d comparisons at %s (%s)"
              % (len(smiles), pairs, a.target, fam), file=sys.stderr)

        ranked = P.rank_ligands(m, smiles, a.target)
        print("%-4s %8s %9s  %-13s %s" % ("RANK", "WIN", "STRENGTH", "BAND IS RIGHT", "COMPOUND"))
        for i, row in enumerate(ranked, 1):
            smi, win, st = row[0], float(row[1]), float(row[2])
            _cut, bacc = band(st, cum)
            label = names.get(smi, smi)
            print("%-4d %8.2f %9.2f  %-13.2f %s" % (i, win, st, bacc, label))
        print("\nHeadline %.2f on %s held-out comparisons, both compounds unseen."
              % (perf["accuracy"], "{:,}".format(perf["held_out_comparisons"])),
              file=sys.stderr)
        return 0

    if a.cmd == "targets":
        rows = [(x, P.family_of(m, x)) for x in P.targets(m)]
        if a.family:
            rows = [x for x in rows if x[1] and a.family.lower() in x[1].lower()]
        for acc, fam in rows:
            print("%-10s %s" % (acc, fam or ""))
        print("\n%d targets" % len(rows), file=sys.stderr)
        return 0

    known = set(P.targets(m))
    asked = [a.a, a.b] if a.cmd == "compare" else a.targets
    unknown = [x for x in asked if x not in known]
    if unknown:
        sys.exit("Not in this model: %s\nTargets are UniProt accessions. "
                 "Run `ffm targets` to list them." % ", ".join(unknown))

    if a.cmd == "compare":
        p = P.compare_targets(m, a.smiles, a.a, a.b)
        s = max(p, 1 - p)
        cut, acc = band(s, cum)
        print("%s (%s)  vs  %s (%s)" % (a.a, P.family_of(m, a.a), a.b, P.family_of(m, a.b)))
        print("  prefers %s at strength %.2f" % (a.a if p >= 0.5 else a.b, s))
        print("  that band is right %.2f of the time" % acc)
        return 0

    ranked = P.rank_targets(m, a.smiles, asked)
    print("%-10s %-38s %8s %9s  %s" % ("TARGET", "FAMILY", "WIN", "STRENGTH", "BAND IS RIGHT"))
    for row in ranked:
        acc_, fam, win, st = row[0], row[1], float(row[2]), float(row[3])
        _cut, bacc = band(st, cum)
        print("%-10s %-38s %8.2f %9.2f  %.2f" % (acc_, (fam or ""), win, st, bacc))
    print("\nHeadline %.2f on %s held-out comparisons."
          % (perf["accuracy"], "{:,}".format(perf["held_out_comparisons"])),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
