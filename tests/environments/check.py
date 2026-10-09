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
the solver is what broke, and only a solve shows that it works.

    PYTHONPATH=src python tests/environments/check.py list
    PYTHONPATH=src python tests/environments/check.py write <dir>
    python tests/environments/check.py run <name> --prefix <conda prefix>

`list` and `write` import CompareM2, so the environment files are the ones
`snakefile.render_envs` generates, never a copy. `run` imports nothing from
CompareM2, so any Python can drive it.
"""

from __future__ import annotations

import argparse
import json
import os
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


EXERCISES = {"basic": fasttree, "carveme": carveme}


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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="check", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="print the environment names as a JSON list")
    w = sub.add_parser("write", help="write <name>.yaml for every environment")
    w.add_argument("directory", type=Path)
    r = sub.add_parser("run", help="exercise one built environment")
    r.add_argument("name", choices=sorted(VERSIONS))
    r.add_argument("--prefix", type=Path, required=True)
    args = p.parse_args(argv)

    if args.command == "list":
        print(json.dumps(sorted(environment_files())))
        return 0
    if args.command == "write":
        args.directory.mkdir(parents=True, exist_ok=True)
        for name, text in environment_files().items():
            (args.directory / f"{name}.yaml").write_text(text)
        return 0
    return run(args.name, args.prefix)


if __name__ == "__main__":
    raise SystemExit(main())
