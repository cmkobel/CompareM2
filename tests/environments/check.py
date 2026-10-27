"""Does each tool environment still build, and does the tool in it still run?

Run weekly by `.github/workflows/environments.yaml`, one environment per job,
because the environments are re-solved by every new install and bioconda moves
underneath them. Two failures made this necessary, and neither was a solve
failure:

- 2026-09-07: the thirteen-tool `main` stopped solving at all (STATUS.md).
- 2026-09-24: `carveme` solved cleanly, installed cleanly, and had no solver.
  scip 10.1.0 arrived under a pyscipopt that links `libscip.so.10.0`, so
  `import pyscipopt` failed and every carveme job died. The published 3.4.0
  was still shipping that on 2026-10-08. A solve-only check passes it.

So every environment is built for real and then exercised. For most tools that
is a version call, which catches a library that will not load or a binary
that is missing, and is all this checks: a tool that crashes only on real
input (bakta 1.8.1 on pyrodigal 3.x) would still pass. `basic` aligns a toy
alignment through FastTree. `carveme` carves a real model through the
pipeline's own `carve_scip.py` and scores it with `biosynthesis.py`, because
the solver is what broke, and only a solve shows that it works. `checkm2`
loads the three models CheckM2 ships and predicts once: its version call
imports TensorFlow but loads no model, and a real run needs the 1.7 GB
database.

The environments are built with conda directly, so none of that goes through
the pipeline. `demo` checks what `comparem2 --demo` produced the way a user
runs it: Snakemake deploying the environments, the generated rules, and the
report. It checks the planted duplicate pair agrees with itself in every tool,
and that a second run finds nothing to do. Panaroo's re-run failure in 3.4.0
was that layer's kind of bug.

`gtdbtk` checks the one database tool CI can run: `--until gtdbtk` against a
40 MB trimmed GTDB package that `gtdb_mock.sh` builds, which classifies the
four test genomes exactly as the real 60.8 GB one does — by ANI screen only,
so prodigal, hmmer and pplacer are never reached.

    PYTHONPATH=src python tests/environments/check.py list
    PYTHONPATH=src python tests/environments/check.py write <dir>
    python tests/environments/check.py run <name> --prefix <conda prefix>
    python tests/environments/check.py demo <output dir> --rerun-log <file>
    python tests/environments/check.py gtdbtk <output dir>

`list` and `write` import CompareM2, so the environment files are the ones
`snakefile.render_envs` generates, never a copy. `run` and `demo` import
nothing from CompareM2, so any Python can drive them.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2] / "src" / "comparem2"

# A version call per tool. Every environment the catalogue declares needs an
# entry, and a unit test holds the two in step.
VERSIONS: dict[str, list[list[str]]] = {
    "basic": [["seqkit", "version"], ["skani", "--version"],
              ["snp-dists", "-v"], ["TreeCluster.py", "--version"],
              ["curl", "--version"], ["tar", "--version"]],
    "perl": [["mlst", "--version"], ["mashtree", "--version"],
             ["panaroo", "--version"]],
    "annotation": [["bakta", "--version"], ["amrfinder", "--version"]],
    "gtdbtk": [["gtdbtk", "--version"]],
    "carveme": [["carve", "--help"]],
    "checkm2": [["checkm2", "--version"]],
}

# The calibration CLAUDE.md records: 29 of 30 de novo for a CarveMe draft of
# E. coli K-12, adenosylcobalamin the one miss. A presolver-mangled model is
# missing hundreds of reactions and lands well below it; re-measure with
# `upstream/panel_calibration.py` before moving this number.
ECOLI_DE_NOVO = 29


def environment_files() -> dict[str, str]:
    """`{name: yaml text}` for the full catalogue, as the pipeline writes them."""
    from comparem2.catalogue import CATALOGUE
    from comparem2.snakefile import render_envs

    return {name.removesuffix(".yaml"): text
            for name, text in render_envs(CATALOGUE, None).items()}


def call(argv: list[str], path: str, **kwargs) -> subprocess.CompletedProcess:
    print("$ " + " ".join(argv), flush=True)
    done = subprocess.run(argv, env={**os.environ, "PATH": path},
                          capture_output=True, text=True, **kwargs)
    tail = (done.stdout + done.stderr).strip().splitlines()[-1:] or ["(no output)"]
    print(f"  exit {done.returncode}: {tail[0][:160]}", flush=True)
    return done


def fasttree(path: str) -> bool:
    toy = ">a\nACGTACGTAC\n>b\nACGTACGTAA\n>c\nACGAACGTAA\n>d\nTCGAACGTAA\n"
    done = call(["FastTree", "-nt", "-gtr", "-quiet"], path, input=toy)
    return done.returncode == 0 and done.stdout.strip().endswith(";")


def carveme(path: str) -> bool:
    """Carve E. coli K-12 from the proteome CarveMe bundles, then score it."""
    where = call(["python", "-c", "import carveme, os; "
                  "print(os.path.dirname(carveme.__file__))"], path)
    if where.returncode:
        return False
    faa = Path(where.stdout.strip()) / "data/benchmark/fasta/Ecoli_K12_MG1655.faa"
    with tempfile.TemporaryDirectory() as tmp:
        # Resolved: `carve_scip.link_input` links the input relatively, and a
        # relative link out of a symlinked /tmp (thylakoid's is /evo/tmp)
        # resolves against the physical parent and dangles. The pipeline never
        # sees this, because its faa and model share a run directory. And
        # check the model exists rather than trusting carve's exit: it exits 0
        # when it cannot open its input.
        tmp = Path(tmp).resolve()
        model, panel = tmp / "ecoli.xml", tmp / "ecoli.tsv"
        if call(["python", str(PACKAGE / "carve_scip.py"), "--faa", str(faa),
                 "--output", str(model)], path).returncode or not model.exists():
            return False
        if call(["python", str(PACKAGE / "biosynthesis.py"), "--model", str(model),
                 "--output", str(panel), "--media", str(tmp / "media.tsv")],
                path).returncode:
            return False
        rows = panel.read_text().splitlines()[1:]
        de_novo = sum(row.split("\t")[3] == "de_novo" for row in rows)
    print(f"  E. coli K-12: {de_novo} of {len(rows)} de novo "
          f"(calibrated at {ECOLI_DE_NOVO})", flush=True)
    return de_novo == ECOLI_DE_NOVO


# What CheckM2's models answer for an all-zero feature vector, measured on
# thylakoid on 2026-10-26 under TensorFlow 2.17.0 and 2.18.0 alike. A zero
# vector is no genome, so the numbers mean nothing biologically; they pin that
# the models a run would use still load and still compute the same function.
CHECKM2_ZERO = "general 6.38 contamination 4.30 specific 2.17"

_CHECKM2_PROBE = """\
import numpy as np
from checkm2.modelProcessing import modelProcessor
m = modelProcessor(1)
x = np.zeros((1, m.general_model_comp.num_feature()))
comp, cont = m.run_prediction_general(x)
nn, _ = m.run_prediction_specific(x, m.specific_model_comp_nn.input_shape[1])
print(f"general {comp[0]:.2f} contamination {cont[0]:.2f} specific {nn[0]:.2f}")
"""


def checkm2(path: str) -> bool:
    """Both LightGBM models, the Keras network and the pickled scaler."""
    done = call(["python", "-c", _CHECKM2_PROBE], path)
    print(f"  zero vector: {done.stdout.strip() or '(nothing)'} "
          f"(measured: {CHECKM2_ZERO})", flush=True)
    return done.returncode == 0 and done.stdout.strip() == CHECKM2_ZERO


EXERCISES = {"basic": fasttree, "carveme": carveme, "checkm2": checkm2}


def run(name: str, prefix: Path) -> int:
    path = f"{prefix / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}"
    failed = [" ".join(argv) for argv in VERSIONS[name]
              if call(argv, path).returncode]
    if name in EXERCISES and not EXERCISES[name](path):
        failed.append(EXERCISES[name].__name__)
    if failed:
        print(f"{name}: FAILED {', '.join(failed)}")
        return 1
    print(f"{name}: ok")
    return 0


# `--demo`'s planted pair: the first plasmid, and the same file again under a
# name with a space, which canonicalisation turns into `116_2_duplicate`.
PAIR = ("116_2", "116_2_duplicate")


def demo(out: Path, rerun_log: Path) -> int:
    """The demo ran all four tools, and the duplicate pair agrees everywhere.

    Values seen on thylakoid on 2026-10-12, from 3.5.0 plus the retry commit:
    the seqkit tables share md5 `e2ec2407…`, skani gives 100.00, and mashtree
    joins the pair with two 0.00000 branches.
    """
    failed = []

    def check(ok: bool, what: str) -> None:
        print(f"  {'ok' if ok else 'FAILED'}: {what}", flush=True)
        if not ok:
            failed.append(what)

    report = (out / "report.html").read_text()
    check("4 of 4 tools produced output" in report, "report: 4 of 4 tools")

    tables = [(out / "samples" / s / "seqkit" / "contigs.tsv").read_bytes() for s in PAIR]
    check(tables[0] == tables[1] and len(tables[0]) > 0, "seqkit: pair identical")

    # skani writes a square matrix: a count, then one row per genome, each
    # named by its path.
    rows = (out / "skani" / "ani.tsv").read_text().splitlines()[1:]
    names = [Path(row.split("\t")[0]).stem for row in rows]
    ani = rows[names.index(PAIR[0])].split("\t")[1:][names.index(PAIR[1])]
    check(float(ani) == 100.0, f"skani: pair at {ani}% ANI")

    newick = (out / "mashtree" / "mashtree.newick").read_text()
    a, b = (re.escape(s) + r":0(\.0+)?" for s in PAIR)
    check(bool(re.search(rf"\(({a},{b}|{b},{a})\)", newick)),
          "mashtree: pair joined at distance 0")

    check("Nothing to be done" in rerun_log.read_text(),
          "re-run: nothing to do")

    print("demo: " + ("FAILED " + ", ".join(failed) if failed else "ok"))
    return 1 if failed else 0


# What GTDB-Tk 2.7.2 says of the four E. faecium test genomes, against the real
# r232 package and the trimmed one `gtdb_mock.sh` builds alike — identical in
# all 20 columns, thylakoid, 2026-10-27. ANI to GCF_029023785.1, the species
# representative.
GTDBTK_SPECIES = "s__Enterococcus_B faecium"
GTDBTK_ANI = {"116_2": 99.18, "116_2_duplicate": 99.18, "E8202": 99.13, "SRR24": 99.03}


def gtdbtk(out: Path) -> int:
    """`--until gtdbtk` against the trimmed package gave the measured answer.

    Only the ANI-screen path: every genome here is classified by skani, so
    GTDB-Tk skips identify and align and the package has no pplacer trees.
    """
    failed = []

    def check(ok: bool, what: str) -> None:
        print(f"  {'ok' if ok else 'FAILED'}: {what}", flush=True)
        if not ok:
            failed.append(what)

    report = (out / "report.html").read_text()
    check("1 of 1 tools produced output" in report, "report: 1 of 1 tools")

    lines = (out / "gtdbtk" / "gtdbtk.summary.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    rows = {r[0]: dict(zip(header, r)) for r in (line.split("\t") for line in lines[1:])}
    check(sorted(rows) == sorted(GTDBTK_ANI), f"summary: {', '.join(sorted(rows))}")
    for genome, ani in GTDBTK_ANI.items():
        row = rows.get(genome, {})
        got = (row.get("classification", "").rsplit(";", 1)[-1],
               row.get("classification_method"), row.get("closest_genome_ani"))
        check(got == (GTDBTK_SPECIES, "ani_screen", f"{ani:.2f}"),
              f"{genome}: {' · '.join(str(g) for g in got)}")

    print("gtdbtk: " + ("FAILED " + ", ".join(failed) if failed else "ok"))
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="check", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="print the environment names as a JSON list")
    w = sub.add_parser("write", help="write <name>.yaml for every environment")
    w.add_argument("directory", type=Path)
    r = sub.add_parser("run", help="exercise one built environment")
    r.add_argument("name", choices=sorted(VERSIONS))
    r.add_argument("--prefix", type=Path, required=True)
    d = sub.add_parser("demo", help="check what `comparem2 --demo` produced")
    d.add_argument("output", type=Path)
    d.add_argument("--rerun-log", type=Path, required=True,
                   help="what a second `--demo` into the same output printed")
    g = sub.add_parser("gtdbtk", help="check `--until gtdbtk` against the trimmed GTDB")
    g.add_argument("output", type=Path)
    args = p.parse_args(argv)

    if args.command == "list":
        print(json.dumps(sorted(environment_files())))
        return 0
    if args.command == "write":
        args.directory.mkdir(parents=True, exist_ok=True)
        for name, text in environment_files().items():
            (args.directory / f"{name}.yaml").write_text(text)
        return 0
    if args.command == "demo":
        return demo(args.output, args.rerun_log)
    if args.command == "gtdbtk":
        return gtdbtk(args.output)
    return run(args.name, args.prefix)


if __name__ == "__main__":
    raise SystemExit(main())
