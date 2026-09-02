"""Run CarveMe's `carve` with SCIP's MILP presolver switched off.

Measured on thylakoid 2026-09-02, one *E. faecium* genome (116_2, 2588
proteins). The MILP is byte-identical in every row — the problem was written
out from each build and the md5s matched — and so is the DIAMOND input:

| SCIP build                             | MILP  | reported  | model              |
| -------------------------------------- | ----- | --------- | ------------------ |
| conda-forge 10.0.3, as shipped         | 601 s | timelimit | 1193 rx / 750 met  |
| the same, `presolving/milp/maxrounds=0`|  43 s | gaplimit  | 1742 rx / 1175 met |
| PyPI wheel 10.0.2 (no PaPILO)          |  10 s | optimal   | 1738 rx / 1176 met |

The difference between the conda-forge build and the wheel is PaPILO, which
only the former links (`printExternalCodeVersions()`). It is not the SCIP
version: conda-forge 10.0.2 carries PaPILO 3.0.0 and is equally slow. It is not
symmetry handling either: `misc/usesymmetry=0` still hit the time limit at
907.8 s.

The shipped run is not only slower, it returns a worse model by CarveMe's own
objective: it drops **253 of 1069 annotated reactions** where the others drop
45. Four configurations of the *one* conda-forge build, on the same written-out
problem, report three different objective values — 913.5 as shipped at the time
limit, 943.5 "optimal" when handed a starting solution, 947.5 "optimal" with
the presolver off — and its dual bound at 300 s is 935.7, below a point the same
build accepts as feasible at 947.5.

**Which of those is the true optimum is not a question this problem answers, and
the fix is not framed as recovering it.** CarveMe couples each flux to its
indicator with bigM=1e3 against eps=1e-3, and every one of those points is
feasible only to tolerance: 5e-11 on the constraints but exactly 1e-6 on
integrality, against SCIP's default `feastol` of 1e-6 — and rounding the
binaries to exact integers makes the remaining LP infeasible in every case. So
"optimal" here describes the solver's tolerances rather than the network. What
is solid is the ordering: with this presolver the run is 12–60x slower and keeps
a quarter fewer of the reactions the annotation supports.

E8202 (3185 proteins) behaves the same way — 908.1 s and 261 annotated reactions
dropped, against 16.7 s and 45 with the presolver off.

Neither CarveMe nor ReFramed exposes solver parameters, which is why this
wrapper exists at all.

**It runs in the tool's environment, under a bare `python`** — the opposite of
`steps.py`, and for the mirror-image reason: that module is CompareM2's and
needs CompareM2's interpreter, this one imports `carveme` and needs the
interpreter that has it. Like `steps.py` it imports nothing from its own
package, because under `--use-conda` its package is not there at all. A test
enforces both halves.

It also hands `carve` an input path inside CarveMe's own output directory.
`carve` derives its DIAMOND output from the *input* path, so pointing it
straight at `bakta/<sample>.faa` writes `bakta/<sample>.tsv` — overwriting
Bakta's feature table with 12-column DIAMOND hits against BiGG. That is what
the 2026-09-02 run left behind, and nothing noticed because Bakta declares only
its GFF3 and its FAA.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# SCIP's PaPILO presolver, disabled. Zero rounds rather than a lower priority:
# the point is that none of its reductions run.
PRESOLVER_OFF = ("presolving/milp/maxrounds", 0)


def link_input(faa: Path, output: Path) -> Path:
    """Link `faa` into `output`'s directory and return the link.

    Relative, so a results directory stays movable, and replaced rather than
    reused so a re-run cannot inherit a link to a genome that has moved.
    """
    outdir = output.parent
    outdir.mkdir(parents=True, exist_ok=True)
    link = outdir / faa.name
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(os.path.relpath(faa, outdir))
    return link


def patch_solver() -> str:
    """Make every SCIP solve skip the MILP presolver. Returns what it did.

    Patched at `solve` rather than at construction because CarveMe sets its own
    `limits/time` and `limits/gap` immediately before solving, so anything set
    earlier is what it overwrites — and this has to be applied after that.

    A missing SCIP backend is reported and left alone. This wrapper exists to
    make a working run faster; it must not be the reason a run fails, and an
    installation solving with Gurobi or CPLEX has nothing here to patch.
    """
    try:
        from reframed.solvers.scip_solver import SCIPSolver
    except Exception as exc:  # ImportError, or pyscipopt failing to load
        return f"SCIP backend unavailable ({exc}); leaving solver parameters alone"

    original = SCIPSolver.solve

    def solve(self, *args, **kwargs):
        self.problem.setParam(*PRESOLVER_OFF)
        return original(self, *args, **kwargs)

    SCIPSolver.solve = solve
    return f"{PRESOLVER_OFF[0]}={PRESOLVER_OFF[1]}"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="carve_scip",
        description="carve, with SCIP's MILP presolver off. See module docstring.")
    p.add_argument("--faa", type=Path, required=True,
                   help="protein FASTA to reconstruct from (Bakta's, normally)")
    p.add_argument("--output", type=Path, required=True, help="SBML model to write")
    args, passthrough = p.parse_known_args(argv)

    print(f"carve_scip: {patch_solver()}", file=sys.stderr)
    link = link_input(args.faa, args.output)

    sys.argv = ["carve", str(link), "--output", str(args.output), *passthrough]
    from carveme.cli.carve import main as carve_main

    return carve_main() or 0


if __name__ == "__main__":
    raise SystemExit(main())
