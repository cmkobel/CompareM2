#!/usr/bin/env python3
"""The biosynthesis verdicts, recomputed with cobrapy and another LP solver.

`biosynthesis.py` reads its verdicts off ReFramed and SCIP. This reads the
same models with cobrapy and solves with HiGHS or GLPK, written from the
module's *documented* definitions rather than from its code:

    de_novo   a drain on M_<id>_c carries flux on M9, every listed compound at
              10 mmol/gDW/h uptake, every other EX_*_e closed to uptake
    upstream  not de novo; carries flux on M9 plus every *other* panel
              compound the model carries
    none      neither
    absent    M_<id>_c is not in the model

and growth, the model's own objective, on each defined reference medium. Only
the panel and the media are taken from the module, as data. Its verdict and
media tables are then compared, verdict for verdict and growth to the four
decimals the table holds. The complete medium is not compared: ReFramed reads
CarveMe's ±1000 bounds as ±∞ (`load_cbmodel(use_infinity=True)`) and cobrapy
keeps them, so the figures there legitimately differ. `iML1515` gives 30.76
against 20.73.

Measured 2026-11-11 on thylakoid: **0 of 1,620 verdicts differ** across 54
models — 12 curated, 42 CarveMe drafts — with HiGHS and again with GLPK, and
every growth figure on the four defined media agrees. `iNF517` cannot be
compared: cobrapy refuses its arginine exchange, fixed at an uptake, once the
medium closes it.

Not in the carveme environment: it needs a Python with cobrapy and a solver,

    python3 -m venv cobra-venv
    cobra-venv/bin/pip install cobra highspy osqp   # osqp: cobrapy's "hybrid"
    cobra-venv/bin/python -I biosynthesis_crosscheck.py \\
        --tables <dir with <model>.tsv and <model>.media.tsv> model.xml ...

`--solver hybrid` is HiGHS; `glpk` is the other one cobrapy ships. Exits 1 on
any disagreement.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path

BIOSYNTHESIS = (Path(__file__).resolve().parent.parent
                / "src" / "comparem2" / "biosynthesis.py")
UPTAKE, MIN_FLUX = 10.0, 1e-6


def load_biosynthesis():
    """The module, for its data. It imports ReFramed only inside functions."""
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Model:
    """A cobrapy model with a drain on every panel compound it carries."""

    def __init__(self, path: Path, panel: list[str]):
        import cobra

        self.cobra = cobra.io.read_sbml_model(str(path))
        self.exchanges = {r.id[3:-2]: r for r in self.cobra.reactions
                          if r.id.startswith("EX_") and r.id.endswith("_e")}
        self.upper = {cid: r.upper_bound for cid, r in self.exchanges.items()}
        self.drains = {}
        for bigg in panel:
            met = f"{bigg}_c"
            if met in self.cobra.metabolites:
                drain = cobra.Reaction(f"XCHECK_DM_{bigg}", lower_bound=0,
                                       upper_bound=1000)
                drain.add_metabolites({self.cobra.metabolites.get_by_id(met): -1})
                self.cobra.add_reactions([drain])
                self.drains[bigg] = drain

    def solve(self, compounds, drain=None) -> tuple[float, str]:
        """Maximise `drain`, or the model's own objective, on `compounds`."""
        with self.cobra:
            for cid, reaction in self.exchanges.items():
                reaction.lower_bound = 0.0
                reaction.upper_bound = self.upper[cid]
            for cid in compounds:
                if cid in self.exchanges:
                    self.exchanges[cid].lower_bound = -UPTAKE
            for other in self.drains.values():
                other.upper_bound = 0.0
            if drain is not None:
                drain.upper_bound = 1000.0
                self.cobra.objective = drain
            value = self.cobra.slim_optimize(error_value=float("nan"))
            status = self.cobra.solver.status
        if status == "unbounded":
            return math.inf, status
        return (value if status == "optimal" else 0.0), status


def verdicts(model: Model, panel: list[str], m9) -> dict[str, str]:
    carried = [b for b in panel if b in model.drains]
    out = {}
    for bigg in panel:
        if bigg not in model.drains:
            out[bigg] = "absent"
        elif model.solve(m9, model.drains[bigg])[0] > MIN_FLUX:
            out[bigg] = "de_novo"
        elif model.solve([*m9, *(b for b in carried if b != bigg)],
                         model.drains[bigg])[0] > MIN_FLUX:
            out[bigg] = "upstream"
        else:
            out[bigg] = "none"
    return out


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def growth_agrees(cell: str, value: float, status: str) -> bool:
    try:
        written = float(cell)
    except ValueError:
        return cell == status  # both a status, e.g. infeasible
    return status == "optimal" and abs(written - max(value, 0.0)) <= 1e-4


def check(bio, path: Path, tables: Path) -> bool:
    panel = [c.bigg for c in bio.PANEL]
    model = Model(path, panel)
    ours = verdicts(model, panel, bio.M9)
    theirs = {r["compound"]: r["verdict"]
              for r in read(tables / f"{path.stem}.tsv")}
    differ = [f"{b}: module {theirs.get(b)} / cobrapy {v}"
              for b, v in ours.items() if theirs.get(b) != v]
    media = {r["medium"]: r["growth"]
             for r in read(tables / f"{path.stem}.media.tsv")}
    for name, compounds, _ in bio.REFERENCE_MEDIA:
        value, status = model.solve(compounds)
        if name in media and not growth_agrees(media[name], value, status):
            differ.append(f"{name} growth: module {media[name]} / "
                          f"cobrapy {value:.4f} {status}")
    counts = {v: sum(x == v for x in ours.values()) for v in sorted(set(ours.values()))}
    print(f"{path.stem:<36} {counts} "
          + ("agrees" if not differ else "DIFFERS: " + "; ".join(differ)),
          flush=True)
    return not differ


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="biosynthesis_crosscheck",
                                description=__doc__.split("\n\n")[0])
    p.add_argument("models", type=Path, nargs="+")
    p.add_argument("--tables", type=Path, required=True,
                   help="where the module's <model>.tsv and .media.tsv are")
    p.add_argument("--solver", default="hybrid", choices=("hybrid", "glpk"))
    args = p.parse_args(argv)

    import cobra

    cobra.Configuration().solver = args.solver
    bio = load_biosynthesis()
    agree = []
    for path in args.models:
        try:
            agree.append(check(bio, path, args.tables))
        except Exception as error:  # a model cobrapy cannot read or solve
            print(f"{path.stem:<36} not compared: {error}", flush=True)
    print(f"{sum(agree)} of {len(agree)} compared models agree "
          f"({len(args.models) - len(agree)} not compared), solver {args.solver}")
    return 0 if all(agree) else 1


if __name__ == "__main__":
    raise SystemExit(main())
