#!/usr/bin/env python3
"""Why does a model report 0 of 32 compounds as de novo?

The biosynthesis panel asks, per compound, whether there is a route to it from
M9 alone. Four of the eight showcase models answer no for all 32 — including
L-alanine and L-glutamate, which are one transamination from central metabolism.
That is either the cascade the panel's own docstring warns about (one missing
link upstream zeroes everything downstream) or a broken model.

Amino acid biosynthesis is almost all transamination from glutamate, and
glutamate comes from 2-oxoglutarate. So if `akg` is unreachable from M9, every
transaminated amino acid is unreachable too and the panel goes to zero for a
single reason rather than 32. This walks central metabolism to find where the
route stops, and then names the reactions that differ between two models.

Self-contained: runs under the tool environment's own python, like
`carve_scip.py` and `biosynthesis.py`, and imports nothing from comparem2.

    python3 probe_akg.py <model.xml> [<model2.xml>]
"""

from __future__ import annotations

import sys

MAX_UPTAKE = 10.0
MIN_FLUX = 1e-6

# The same M9 as the panel: salts, glucose, ammonium, phosphate, sulfate, O2.
M9 = ("ca2", "cl", "co2", "cobalt2", "cu2", "fe2", "fe3", "glc__D", "h", "h2o",
      "k", "mg2", "mn2", "mobd", "na1", "nh4", "ni2", "o2", "pi", "so4", "zn4",
      "zn2")

# Glycolysis down to pyruvate, then the routes that would reach 2-oxoglutarate,
# then the amino acids that depend on it. Ordered so the first `no` is the
# bottleneck.
PROBES = [
    ("glc__D", "glucose (uptake control)"),
    ("g6p", "glucose-6-phosphate"),
    ("3pg", "3-phosphoglycerate"),
    ("pep", "phosphoenolpyruvate"),
    ("pyr", "pyruvate"),
    ("accoa", "acetyl-CoA"),
    ("oaa", "oxaloacetate"),
    ("cit", "citrate"),
    ("icit", "isocitrate"),
    ("akg", "2-oxoglutarate  <- the hinge"),
    ("glu__L", "L-glutamate"),
    ("gln__L", "L-glutamine"),
    ("ala__L", "L-alanine"),
    ("asp__L", "L-aspartate"),
    ("ser__L", "L-serine"),
]


def load(path: str):
    from reframed import load_cbmodel

    return load_cbmodel(path, flavor="bigg")


def exchanges(model) -> dict[str, str]:
    return {rid[len("R_EX_"):-len("_e")]: rid for rid in model.reactions
            if rid.startswith("R_EX_") and rid.endswith("_e")}


def producible(path: str) -> dict[str, float]:
    """Max flux through a drain on each probe metabolite, on M9."""
    from reframed import FBA
    from reframed.solvers import solver_instance

    model = load(path)
    ex = exchanges(model)

    # One drain per probe present, added before the solver is built.
    drains = {}
    for bigg, _ in PROBES:
        if f"M_{bigg}_c" not in model.metabolites:
            continue
        rid = f"R_DM_probe_{bigg}"
        model.add_reaction_from_str(f"{rid}: M_{bigg}_c --> ")
        drains[bigg] = rid

    solver = solver_instance(model)
    upper = {rid: model.reactions[rid].ub for rid in ex.values()}

    out = {}
    for bigg, rid in drains.items():
        # Exactly M9 available; every other exchange shut; every probe drain
        # shut but this one, so no drain can act as a free sink for another.
        cons = {r: (0.0, upper[r]) for r in ex.values()}
        cons.update({r: (0.0, 0.0) for r in drains.values()})
        for cid in M9:
            if cid in ex:
                cons[ex[cid]] = (-MAX_UPTAKE, upper[ex[cid]])
        cons[rid] = (0.0, 1000.0)
        sol = FBA(model, objective={rid: 1}, constraints=cons, solver=solver)
        status = str(getattr(sol.status, "value", sol.status))
        out[bigg] = (sol.fobj or 0.0) if status == "Optimal" else 0.0
    return out


def producers(path: str, bigg: str) -> set[str]:
    """Reactions that can carry `bigg` on the product side."""
    model = load(path)
    met = f"M_{bigg}_c"
    found = set()
    for rid, reaction in model.reactions.items():
        coeff = reaction.stoichiometry.get(met)
        if coeff is None:
            continue
        # Either direction counts if the reaction is reversible.
        if coeff > 0 or reaction.lb < 0:
            found.add(rid)
    return found


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        raise SystemExit(__doc__)
    paths = argv[1:]

    results = {}
    for path in paths:
        name = path.rsplit("/", 1)[-1]
        results[name] = producible(path)
        print(f"loaded {name}")

    print()
    width = max(len(label) for _, label in PROBES) + 2
    header = "metabolite".ljust(12) + "route from M9".ljust(width)
    header += "".join(n.ljust(22) for n in results)
    print(header)
    for bigg, label in PROBES:
        row = bigg.ljust(12) + label.ljust(width)
        for name in results:
            value = results[name].get(bigg)
            if value is None:
                cell = "absent"
            else:
                cell = f"{'YES' if value > MIN_FLUX else 'no':4s} {value:8.3f}"
            row += cell.ljust(22)
        print(row)

    if len(paths) == 2:
        for bigg in ("akg", "glu__L"):
            a, b = (producers(p, bigg) for p in paths)
            na, nb = (p.rsplit("/", 1)[-1] for p in paths)
            print(f"\nreactions touching {bigg}: {len(a)} in {na}, "
                  f"{len(b)} in {nb}")
            only_a, only_b = sorted(a - b), sorted(b - a)
            if only_a:
                print(f"  only in {na}: {', '.join(only_a)}")
            if only_b:
                print(f"  only in {nb}: {', '.join(only_b)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
