#!/usr/bin/env python3
"""Negative controls for the biosynthesis probe: what a model must *not* make.

The panel is a list of things a model can build. A model can also build things
it has no right to — a mass-imbalanced reaction, an open sink, an exchange the
medium does not govern — and every such route shows up as a *de novo* call that
looks exactly like a real one. This asks each model for the panel on media that
lack an element the compounds need:

    empty     nothing at all          nothing may be made
    M9-C      M9 without glucose      nothing: every panel compound has carbon
    M9-N      M9 without ammonium     only menaquinone-8, the one without nitrogen
    M9-S      M9 without sulfate      not cysteine, methionine, thiamine diphosphate
    M9-P      M9 without phosphate    not NAD, ThDP, PLP, adenosylcobalamin
    M9-Fe     M9 without iron         not protoheme, siroheme
    M9-Co     M9 without cobalt       not adenosylcobalamin

A violation is a compound made on a medium whose missing element it contains.
It also reports what the medium does not govern: one-sided reactions other than
`R_EX_*_e` that could supply, forced fluxes, and exchanges not written
`x_e <=>`.

Measured 2026-11-11 on thylakoid, 42 CarveMe drafts — 17 from pipeline runs,
the six 1.5.0 benchmark drafts, 19 carved that day: **no violation, no
supplier, no forced flux, every exchange oriented.** Four of 13 curated BiGG
models leak. `iIT341` makes cysteine with no sulfate through its reversible
`SK_ahcys_c` sink. The other three keep their leak with every sink shut,
because it is a reaction that does not balance. `iYS854` makes 22 panel
compounds with no nitrogen; traced for alanine, the nitrogen comes from `SRK`,
which makes phosphoserine from 3-phosphohydroxypyruvate with no amino donor.
`iCN900` makes NAD with no phosphate, through a `GLNSCP` written with one ATP.

Self-contained, like `panel_calibration.py`: run it under the carveme
environment's python, which has ReFramed. It imports `biosynthesis.py` by path.

    python3 biosynthesis_controls.py model.xml [model.xml ...]

Exits 1 if any model has a violation.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

BIOSYNTHESIS = (Path(__file__).resolve().parent.parent
                / "src" / "comparem2" / "biosynthesis.py")

# Compounds a medium without the element must not yield. `None` is all of them.
FORBIDDEN = {
    "empty": None,
    "M9-C": None,
    "M9-N": "all but mqn8",
    "M9-S": {"cys__L", "met__L", "thmpp"},
    "M9-P": {"nad", "thmpp", "pydx5p", "adocbl"},
    "M9-Fe": {"pheme", "sheme"},
    "M9-Co": {"adocbl"},
}


def load_biosynthesis():
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def control_media(bio) -> dict[str, tuple[str, ...]]:
    def without(*drop):
        return tuple(c for c in bio.M9 if c not in drop)

    return {"empty": (), "M9-C": without("glc__D"), "M9-N": without("nh4"),
            "M9-S": without("so4"), "M9-P": without("pi"),
            "M9-Fe": without("fe2", "fe3"), "M9-Co": without("cobalt2")}


def forbidden(label: str, bigg: str) -> bool:
    rule = FORBIDDEN[label]
    if rule is None:
        return True
    if rule == "all but mqn8":
        return bigg != "mqn8"
    return bigg in rule


def suppliers(model, exchanges) -> list[str]:
    """One-sided reactions other than the exchanges that can supply a metabolite.

    `x <=>` supplies when it runs backwards, `<=> x` when it runs forwards.
    """
    out = []
    for rid, reaction in model.reactions.items():
        coefficients = list(reaction.stoichiometry.values())
        if rid in exchanges or not coefficients:
            continue
        lb, ub = reaction.lb, reaction.ub  # None is unbounded
        consumes = all(c < 0 for c in coefficients)
        produces = all(c > 0 for c in coefficients)
        if (consumes and (lb is None or lb < 0)) or (produces and (ub is None or ub > 0)):
            out.append(rid)
    return sorted(out)


def audit(bio, model) -> dict[str, list]:
    """What the medium does not govern, read off the model before any drain."""
    exchanges = set(bio._exchanges(model).values())
    misoriented = [rid for rid in exchanges
                   if dict(model.reactions[rid].stoichiometry)
                   != {f"M_{rid[len('R_EX_'):-len('_e')]}_e": -1}]
    forced = [rid for rid, r in model.reactions.items()
              if (r.lb is not None and r.lb > 0) or (r.ub is not None and r.ub < 0)]
    return {"suppliers": suppliers(model, exchanges),
            "forced": sorted(forced), "misoriented": sorted(misoriented)}


def run(bio, path: Path) -> bool:
    from reframed import load_cbmodel

    model = load_cbmodel(str(path), flavor="bigg")
    found = audit(bio, model)
    probe = bio._Probe(model)
    violations = []
    for label, compounds in control_media(bio).items():
        for compound in probe.present:
            if (forbidden(label, compound.bigg)
                    and probe.reaches(bio.cytoplasmic(compound.bigg), compounds)):
                violations.append(f"{label}:{compound.bigg}")
    print(f"{path.stem:<36} suppliers={len(found['suppliers'])} "
          f"forced={len(found['forced'])} misoriented={len(found['misoriented'])} "
          f"statuses={probe.statuses} "
          f"violations={' '.join(violations) or 'none'}", flush=True)
    for kind in ("suppliers", "forced", "misoriented"):
        if found[kind]:
            print(f"    {kind}: {' '.join(found[kind][:8])}"
                  + (" …" if len(found[kind]) > 8 else ""), flush=True)
    return not violations


def main(argv: list[str] | None = None) -> int:
    paths = [Path(a) for a in (argv if argv is not None else sys.argv[1:])]
    if not paths:
        print(__doc__, file=sys.stderr)
        return 2
    bio = load_biosynthesis()
    clean = [run(bio, path) for path in paths]
    print(f"{sum(clean)} of {len(clean)} models clean")
    return 0 if all(clean) else 1


if __name__ == "__main__":
    raise SystemExit(main())
