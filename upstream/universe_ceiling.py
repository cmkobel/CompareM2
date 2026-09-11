#!/usr/bin/env python3
"""Can any CarveMe draft ever answer `de novo` for this compound?

Every draft is a subnetwork of the BiGG universe CarveMe carves from, so a
compound the universe cannot reach from M9 cannot be reached by any draft. That
makes this the ceiling on the biosynthesis panel, and the check that separates
"the database has no route" from "carving did not keep the route".

**The answer on 2026-09-10 was that there is no ceiling**: all 32 compounds on
the panel as it then stood are reachable, so every probe target is a
representation the database can actually produce and no verdict is an artefact
of a badly chosen metabolite id. Which also settled the other question — `btn`
and `q8` came out `none` in 12 of 12 drafts across 8 species, and since the
routes are in the universe, that is carving and not the database. Both were
dropped from the panel the same day. Re-run this when CarveMe's universe moves.

Two things the universe needs before it can be asked.

- **It ships without exchange reactions.** `carve` adds them, so the M9 ones
  are added here.
- **Every one of its 25,348 reactions is bounded at ±inf**, which makes the LP
  unbounded rather than optimal. That mattered: before `biosynthesis.py` learned
  to read `Unbounded` as producible, this script reported carbon *and* nitrogen
  unreachable and 29 of 32 compounds `none` — every one of them an unbounded
  solve misread as zero flux. The bounds are capped here anyway, because an
  unbounded network answers a weaker question than a bounded one, and the status
  tally is printed either way.

Self-contained, like `carve_scip.py` and `biosynthesis.py` — runs under the tool
environment's own python and imports nothing from `comparem2`.

    python3 universe_ceiling.py
    python3 universe_ceiling.py --universe /path/to/bigg_universe.xml
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIOSYNTHESIS = HERE.parent / "src" / "comparem2" / "biosynthesis.py"

# The universe ships unbounded; anything finite and large enough not to bind
# will do, and 1000 is what every BiGG model uses for "no limit".
BOUND = 1000.0


def load_biosynthesis():
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bundled_universe() -> Path:
    """Where CarveMe keeps the universe in the installed package."""
    import carveme

    return (Path(carveme.__file__).parent / "data" / "generated"
            / "bigg_universe.xml.gz")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="universe_ceiling", description=__doc__)
    p.add_argument("--universe", type=Path, default=None,
                   help="SBML universe (default: the one CarveMe installed)")
    p.add_argument("--bound", type=float, default=BOUND,
                   help=f"cap every reaction at ±this (default {BOUND})")
    args = p.parse_args(argv)

    bio = load_biosynthesis()
    from reframed import load_cbmodel

    path = args.universe or bundled_universe()
    with tempfile.TemporaryDirectory() as tmp:
        if path.suffix == ".gz":
            plain = Path(tmp) / path.stem
            with gzip.open(path) as src, plain.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            path = plain
        model = load_cbmodel(str(path), flavor="bigg")

    print(f"universe: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites", file=sys.stderr)
    for reaction in model.reactions.values():
        reaction.lb, reaction.ub = -args.bound, args.bound
    for compound in bio.M9:
        model.add_reaction_from_str(f"R_EX_{compound}_e: M_{compound}_e <-> ")

    probe = bio._Probe(model)
    print(f"panel compounds present: {len(probe.present)}/{len(bio.PANEL)}")
    print(f"unreachable source elements on M9: "
          f"{bio.unreachable_sources(probe, bio.M9) or 'none'}")

    by: dict[str, list[str]] = {}
    for bigg, _, _, verdict in bio.verdicts(probe):
        by.setdefault(verdict, []).append(bigg)
    for verdict in (bio.DE_NOVO, bio.UPSTREAM, bio.NO_ROUTE, bio.ABSENT):
        got = by.get(verdict, [])
        print(f"  {verdict:>8} {len(got):>2}: {' '.join(got)}")
    print(f"statuses: {probe.statuses}")

    impossible = by.get(bio.NO_ROUTE, []) + by.get(bio.ABSENT, [])
    print("\nout of reach for every CarveMe draft: "
          + (" ".join(impossible) if impossible else "none — there is no ceiling"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
