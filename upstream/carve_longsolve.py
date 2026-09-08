#!/usr/bin/env python3
"""Re-carve one genome with SCIP's time limit raised, to separate two causes.

The showcase run produced two models of *the same strain* — `Spn_D39` and its
unencapsulated derivative `Spn_R6`, 65 core SNPs apart — that disagree on 20 of
the 32 biosynthesis panel compounds, D39 `de_novo` where R6 says `upstream`.
D39's solve stopped at CarveMe's 600 s ceiling with a 2.6% gap; R6's returned
`optimal` in 15.4 s. So there are two candidate explanations and they call for
different responses:

  truncation  D39 is simply unfinished. Given time it converges on R6's answer,
              and the fix is solver time — the report is slow, not wrong.
  degeneracy  CarveMe's near-optimal region holds networks with different
              biosynthetic closure and the search path picks one. Then no time
              limit fixes it and the panel cannot be read per genome at all.

`Spn_P1031` is the reason this is not obvious: it met CarveMe's own gap
criterion (7.6e-4, under `limits/gap=0.001`) and still answers 18 de novo like
the truncated models, against 0 for the four other converged ones.

Same patch point as `carve_scip.py` — `SCIPSolver.solve`, because CarveMe sets
`limits/time` and `limits/gap` immediately before solving and anything set
earlier is overwritten. This is a one-off measurement, not a pipeline change;
nothing in `src/` is touched.

    python3 carve_longsolve.py --faa <bakta.faa> --output <dir>/<name>.xml \
        --seconds 7200 [--gap 1e-4]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PRESOLVER_OFF = ("presolving/milp/maxrounds", 0)

SOLVE_FIELDS = (
    ("status", "getStatus"),
    ("gap", "getGap"),
    ("primal", "getPrimalbound"),
    ("dual", "getDualbound"),
    ("seconds", "getSolvingTime"),
    ("solutions", "getNSols"),
)


def describe(problem) -> str:
    parts = []
    for label, getter in SOLVE_FIELDS:
        try:
            value = getattr(problem, getter)()
            if isinstance(value, float):
                value = f"{value:.6g}"
        except Exception as exc:
            value = f"unavailable({type(exc).__name__})"
        parts.append(f"{label}={value}")
    return " ".join(parts)


def patch(seconds: float, gap: float | None) -> None:
    from reframed.solvers.scip_solver import SCIPSolver

    original = SCIPSolver.solve

    def solve(self, *args, **kwargs):
        try:
            self.problem.setParam(*PRESOLVER_OFF)
        except KeyError:
            print("longsolve: no PaPILO in this build", file=sys.stderr)
        # After CarveMe's own settings, which is the whole point of patching
        # here rather than at construction.
        self.problem.setParam("limits/time", seconds)
        if gap is not None:
            self.problem.setParam("limits/gap", gap)
        solution = original(self, *args, **kwargs)
        print(f"longsolve: solve {describe(self.problem)}", file=sys.stderr)
        return solution

    SCIPSolver.solve = solve


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--faa", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seconds", type=float, default=7200.0)
    p.add_argument("--gap", type=float, default=None)
    args, passthrough = p.parse_known_args(argv)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # carve derives its DIAMOND output path from its *input* path, so the FAA is
    # linked next to the model rather than read in place — otherwise the hits
    # land on top of Bakta's feature table. Same reason as carve_scip.py.
    link = args.output.parent / args.faa.name
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(args.faa.resolve())

    print(f"longsolve: limits/time={args.seconds} limits/gap={args.gap}",
          file=sys.stderr)
    patch(args.seconds, args.gap)

    sys.argv = ["carve", str(link), "--output", str(args.output), *passthrough]
    from carveme.cli.carve import main as carve_main

    return carve_main() or 0


if __name__ == "__main__":
    raise SystemExit(main())
