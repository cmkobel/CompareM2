#!/usr/bin/env python3
"""What the biosynthesis panel answers for genomes whose answer is known.

The panel's calibration used to be one curated model, `iML1515` at 31 of 32
de novo. That is the wrong reference class: this pipeline never produces a
curated model, and on 2026-09-10 a CarveMe draft of the *same organism* came
out at 29, missing biotin and ubiquinone-8 — both of which *E. coli* makes.
Dropping those two took the panel to 30 and the curated model and the drafts
then agreed at **29 of 30, adenosylcobalamin the only miss in either**.

This re-measures that. It is the check to re-run when CarveMe's version, its
universe or its scoring changes, because the agreement is what says how much
weight the *de novo* column carries.

**The inputs ship with CarveMe.** `carveme/data/benchmark/fasta/` holds the
proteomes of the five organisms the CarveMe paper benchmarks against plus
*M. genitalium*, so nothing has to be downloaded and nothing has to be
annotated. Carving *E. coli* K-12 took 7.5 s of solver time on thylakoid and
*B. subtilis* 168 took 23.7 s, both to a certified optimum;
*M. genitalium* took 274.9 s and stopped at a 6.9e-5 gap.

Measured 2026-09-10, CarveMe 1.6.6 and ReFramed 1.6.0, panel of 30:

    iML1515             29 de novo   adocbl none    (curated, for comparison)
    Ecoli_K12_MG1655    29 de novo   adocbl absent
    Bsubtilis_168       29 de novo   adocbl absent
    M_genitalium_G37     0 de novo   3 upstream, 24 none, 3 absent

The *M. genitalium* row is the negative control and the reason it is in the
set: a genome-reduced obligate parasite must not come out prototrophic.

Self-contained, like `carve_scip.py` and `biosynthesis.py` — runs under the
tool environment's own python. It imports `biosynthesis` by path rather than
from `comparem2`, because under `--use-conda` the package is not installed in
the environment that has ReFramed.

**Run it with that environment activated, not just with its interpreter**:
`carve` shells out to DIAMOND, so the environment's `bin` has to be on `PATH`
or the carve step fails with "Unable to run diamond". Models already in the
workdir are reused, so scoring alone needs no DIAMOND at all.

    conda activate <the carveme env>
    python3 panel_calibration.py --workdir /tmp/cal
    python3 panel_calibration.py --workdir /tmp/cal --model iML1515.xml
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIOSYNTHESIS = HERE.parent / "src" / "comparem2" / "biosynthesis.py"
CARVE_SCIP = HERE.parent / "src" / "comparem2" / "carve_scip.py"

# The prototroph controls, then the negative control. Named as the files are.
GENOMES = ("Ecoli_K12_MG1655", "Bsubtilis_168", "M_genitalium_G37")


def load_biosynthesis():
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bundled_fasta() -> Path:
    """Where CarveMe keeps the benchmark proteomes in the installed package."""
    import carveme

    return Path(carveme.__file__).parent / "data" / "benchmark" / "fasta"


def carve(faa: Path, output: Path) -> None:
    """Through `carve_scip.py`, so this measures the pipeline's own path."""
    if output.exists():
        print(f"  {output.name} exists, reusing", file=sys.stderr)
        return
    subprocess.run([sys.executable, str(CARVE_SCIP),
                    "--faa", str(faa), "--output", str(output)], check=True)


def report(bio, path: Path) -> tuple[str, dict[str, int], list[str]]:
    from reframed import load_cbmodel

    model = load_cbmodel(str(path), flavor="bigg")
    probe = bio._Probe(model)
    counts: dict[str, int] = {}
    missed = []
    for compound, _, _, verdict in bio.verdicts(probe):
        counts[verdict] = counts.get(verdict, 0) + 1
        if verdict != bio.DE_NOVO:
            missed.append(f"{compound}:{verdict}")
    return path.stem, counts, missed


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="panel_calibration", description=__doc__)
    p.add_argument("--workdir", type=Path, required=True,
                   help="where to carve into; existing models are reused")
    p.add_argument("--model", type=Path, action="append", default=[],
                   help="an extra model to score, e.g. a curated one. Repeatable")
    p.add_argument("--genome", action="append", default=[],
                   help=f"override the bundled set (default: {', '.join(GENOMES)})")
    args = p.parse_args(argv)

    bio = load_biosynthesis()
    args.workdir.mkdir(parents=True, exist_ok=True)

    models = list(args.model)
    fasta = bundled_fasta()
    for name in (args.genome or GENOMES):
        faa = fasta / f"{name}.faa"
        if not faa.exists():
            print(f"no bundled proteome {faa}", file=sys.stderr)
            return 1
        # Into the workdir, and not next to the proteome: `carve_scip.py` links
        # its input into the output's directory.
        output = args.workdir / f"{name}.xml"
        carve(faa, output)
        models.append(output)

    print(f"{'model':<22} {'panel':>5} {'de novo':>8}  not de novo")
    for path in models:
        name, counts, missed = report(bio, Path(path))
        print(f"{name:<22} {len(bio.PANEL):>5} "
              f"{counts.get(bio.DE_NOVO, 0):>8}  {' '.join(missed) or '—'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
