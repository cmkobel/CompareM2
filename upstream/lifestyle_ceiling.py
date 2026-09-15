#!/usr/bin/env python3
"""Is a lifestyle out of reach for every CarveMe draft, whatever the scoring?

A draft is a subnetwork of the universe it was carved from, so a reaction the
universe does not contain cannot appear in any draft of any genome. That makes
this the ceiling on a lifestyle table, and the one experiment that can kill the
idea outright rather than merely disappoint it — `upstream/lifestyle.md` lists
it first for exactly that reason.

It settles the inference the whole plan rested on and which had never been
measured: *that reactions present in a BiGG model survive into the universe
CarveMe carves from.* The 2026-09-15 entry in `DECISIONS.md` called it
"reasonable and not a measurement". Measured, **it is false for precisely the
reactions the idea needs**.

## The answer, measured 2026-09-15 on thylakoid

carveme 1.6.6, reframed 1.6.0.

    universe_bacteria       5,532 reactions
    universe_gramneg        5,571
    universe_grampos        5,664
    universe_cyanobacteria  5,680
    universe_archaea        5,724
    bigg_universe          25,348
    union of the five       5,948   — 20,070 of BiGG's reactions are in none of them

**First correction: the carving universe is 5,532 reactions, not 25,348.** That
figure is `bigg_universe.xml.gz`, an input to CarveMe's build rather than a
thing `carve` ever loads; `config.cfg` points `default_universe` at
`universe_bacteria.xml.gz`. Several places in this repo used the larger number
for the thing being carved from, and arm three's cost estimate was four times
too high as a result.

Every one of the nineteen marker reactions below is in `bigg_universe`. What
they are missing from is the universe that is actually used.

| lifestyle | marker | bact | arch | cyano |
| --- | --- | :-: | :-: | :-: |
| photoautotrophy | `RBPC` | no | yes | **yes** |
| | `RBCh` `PSI` `PSII` `PSI_2` | no | no | **yes** |
| Wood–Ljungdahl | `FTHFLi` `MTHFC` `MTHFD` | yes | yes | yes |
| | `CODH_ACS` `CODH4` `MTHFR5` `RNF` `HYDFDN2r` `FDH7` | **no** | **no** | **no** |
| methanogenesis | `MCR` `HDR` `FMFD_b` `CODHr` `CODH2r` | **no** | **no** | **no** |
| the acceptor axis | `FRD*` `NO3R*` `DMSOR*` `TMAOR*` `SULR*` `FE3Ri` | yes | yes | yes |

**Second correction, and it reverses a decision taken the same day.** The
archaeal universe does **not** restore methanogenesis. `--set
carveme--universe=archaea` was adopted on 2026-09-15 as the cheap fix for
`MCR` and `HDR` living in an archaeal model — and `universe_archaea` does not
contain them either. Checked at the chemistry level so it cannot be a naming
difference: **coenzyme M, coenzyme B and methanofuran are not metabolites of any
of the five universes.** Methanogenesis is not representable, not merely
unscored. `M_ch4_c` exists in `universe_archaea` with nothing that makes it.

What the archaeal universe does buy is a better fit generally — `iAF692`
overlaps it at 542 of 690 reactions against 431 for bacteria — so the flag is
still right for archaea. It just does not buy the column it was adopted for.

**Third: photoautotrophy is reachable, from `universe_cyanobacteria` only.** All
five of `iJN678`'s photosynthesis reactions are there and four of the five are
in no other universe. That was not known and no decision had been taken on it.

So the ceiling is: **autotrophy, lithotrophy and methanogenesis cannot be
reached by any CarveMe draft of any genome**, which is a database limit and the
plan's own test — "if the universe cannot reach CH4 from H2 + CO2, that is a
database limit and no scoring rescues it". The acceptor axis is fully supported.
Those are the two halves the plan expected to be the other way round.

Self-contained, like `carve_scip.py` and `universe_ceiling.py` — the tool
environment's own python, nothing imported from `comparem2`.

    conda activate <the carveme env>
    python3 lifestyle_ceiling.py --workdir /tmp/ceiling
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

BIGG_URL = "https://bigg.ucsd.edu/static/models/{model}.xml.gz"

UNIVERSES = ("bacteria", "gramneg", "grampos", "cyanobacteria", "archaea")

# The marker reactions, grouped by the curated model they were verified in on
# 2026-09-15 (`DECISIONS.md`). **None is from memory** — each was read out of the
# downloaded SBML, which is why this list and not a longer, guessed one.
MARKERS = (
    ("iJN678", "photoautotrophy",
     ("R_RBPC", "R_RBCh", "R_PSI", "R_PSII", "R_PSI_2")),
    ("iHN637", "Wood-Ljungdahl, Rnf, bifurcating hydrogenase",
     ("R_CODH_ACS", "R_CODH4", "R_FTHFLi", "R_MTHFC", "R_MTHFD", "R_MTHFR5",
      "R_RNF", "R_HYDFDN2r", "R_FDH7")),
    ("iAF692", "methanogenesis",
     ("R_MCR", "R_HDR", "R_FMFD_b", "R_CODHr", "R_CODH2r")),
)

# Asked as chemistry rather than as reaction ids, because a reaction absent
# under one name could be present under another — a metabolite that is not in
# the universe at all cannot be. These are the cofactors methanogenesis runs on.
COFACTORS = ("com", "cob", "mfr_b", "ch4", "h2", "co")

# The acceptor axis, matched by prefix rather than by exact id: the universes
# carry periplasmic and cytosolic variants under several suffixes, and what
# matters is whether *any* of them is there.
ACCEPTORS = (("fumarate", "FRD"), ("nitrate", "NO3R"), ("DMSO", "DMSOR"),
             ("TMAO", "TMAOR"), ("sulfite/sulfate", "SULR"),
             ("Fe(III)", "FE3R"))


def load(path: Path):
    """Read an SBML model, gunzipping into a scratch file if it needs it."""
    from reframed import load_cbmodel

    with tempfile.TemporaryDirectory() as tmp:
        plain = Path(tmp) / "model.xml"
        if path.suffix == ".gz":
            with gzip.open(path) as src, plain.open("wb") as dst:
                shutil.copyfileobj(src, dst)
        else:
            shutil.copy(path, plain)
        return load_cbmodel(str(plain), flavor="bigg")


def generated() -> Path:
    """Where CarveMe keeps the universes in the installed package."""
    import carveme

    return Path(carveme.__file__).parent / "data" / "generated"


def curated(model_id: str, workdir: Path) -> Path:
    """One BiGG model, downloaded once."""
    path = workdir / f"{model_id}.xml"
    if not path.exists() or not path.stat().st_size:
        print(f"  fetching {model_id}", file=sys.stderr)
        body = urllib.request.urlopen(BIGG_URL.format(model=model_id)).read()
        path.write_bytes(gzip.decompress(body))
    return path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="lifestyle_ceiling",
                                description=__doc__.split("\n\n")[0])
    p.add_argument("--workdir", type=Path, required=True,
                   help="where to download curated models; existing files reused")
    args = p.parse_args(argv)
    args.workdir.mkdir(parents=True, exist_ok=True)

    folder = generated()
    universes = {name: load(folder / f"universe_{name}.xml.gz")
                 for name in UNIVERSES}
    big = load(folder / "bigg_universe.xml.gz")

    print("universe sizes")
    for name in UNIVERSES:
        print(f"  universe_{name:<15} {len(universes[name].reactions):>6} reactions")
    print(f"  {'bigg_universe':<24} {len(big.reactions):>6} reactions")
    union = set().union(*(set(u.reactions) for u in universes.values()))
    print(f"  {'union of the five':<24} {len(union):>6}")
    print(f"  {'in BiGG, in no universe':<24} "
          f"{len(set(big.reactions) - union):>6}\n")

    names = list(UNIVERSES)
    head = "  ".join(f"{n[:5]:>5}" for n in names)
    print(f"{'marker':<14} {'model':>5}  {head}  {'bigg':>5}")
    for model_id, label, reactions in MARKERS:
        model = load(curated(model_id, args.workdir))
        overlap = {n: len(set(model.reactions) & set(u.reactions))
                   for n, u in universes.items()}
        print(f"-- {model_id} — {label}: {len(model.reactions)} reactions, "
              + ", ".join(f"{n} {overlap[n]}" for n in names))
        for rid in reactions:
            cells = "  ".join(
                f"{('yes' if rid in universes[n].reactions else '-'):>5}"
                for n in names)
            print(f"{rid:<14} {('yes' if rid in model.reactions else '-'):>5}  "
                  f"{cells}  {('yes' if rid in big.reactions else '-'):>5}")

    print(f"\n{'cofactor':<14} {'':>5}  {head}  {'bigg':>5}")
    for met in COFACTORS:
        mid = f"M_{met}_c"
        cells = "  ".join(
            f"{('yes' if mid in universes[n].metabolites else '-'):>5}"
            for n in names)
        print(f"{met:<14} {'':>5}  {cells}  "
              f"{('yes' if mid in big.metabolites else '-'):>5}")

    print(f"\n{'acceptor':<14} {'':>5}  {head}")
    for label, prefix in ACCEPTORS:
        cells = "  ".join(
            f"{sum(1 for r in universes[n].reactions if r[2:].startswith(prefix)):>5}"
            for n in names)
        print(f"{label:<14} {'':>5}  {cells}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
