#!/usr/bin/env python3
"""Which reaction lets one model fix nitrogen into amino acids, and why is it there?

`probe_akg.py` showed both showcase models unable to reach 2-oxoglutarate or
glutamate from M9, yet `Spn_D39` reaches alanine, aspartate and serine and its
near-identical sibling `Spn_R6` does not. Those are normally transaminations
from glutamate, so D39 has a glutamate-independent entry for ammonium.

This finds it: maximise a drain on aspartate under M9, list the reactions
carrying flux that consume ammonium, and report whether each has a gene
association in the SBML. That last column is the point. A reaction placed by a
DIAMOND hit to an annotated gene is evidence; a reaction with no gene
association is CarveMe's gap-filling, and a gap-fill that invents an ammonium
entry point is what would make 18 de novo calls spurious.

    python3 find_entry.py <model.xml> [<metabolite>]
"""

from __future__ import annotations

import sys

MAX_UPTAKE = 10.0
MIN_FLUX = 1e-6

M9 = ("ca2", "cl", "co2", "cobalt2", "cu2", "fe2", "fe3", "glc__D", "h", "h2o",
      "k", "mg2", "mn2", "mobd", "na1", "nh4", "ni2", "o2", "pi", "so4", "zn4",
      "zn2")


def main(argv: list[str]) -> int:
    from reframed import FBA, load_cbmodel
    from reframed.solvers import solver_instance

    if not 2 <= len(argv) <= 3:
        raise SystemExit(__doc__)
    path = argv[1]
    target = argv[2] if len(argv) == 3 else "asp__L"

    model = load_cbmodel(path, flavor="bigg")
    ex = {rid[len("R_EX_"):-len("_e")]: rid for rid in model.reactions
          if rid.startswith("R_EX_") and rid.endswith("_e")}

    drain = f"R_DM_probe_{target}"
    model.add_reaction_from_str(f"{drain}: M_{target}_c --> ")
    solver = solver_instance(model)

    upper = {rid: model.reactions[rid].ub for rid in ex.values()}
    cons = {rid: (0.0, upper[rid]) for rid in ex.values()}
    for cid in M9:
        if cid in ex:
            cons[ex[cid]] = (-MAX_UPTAKE, upper[ex[cid]])
    cons[drain] = (0.0, 1000.0)

    sol = FBA(model, objective={drain: 1}, constraints=cons, solver=solver)
    status = str(getattr(sol.status, "value", sol.status))
    print(f"{path.rsplit('/', 1)[-1]}: max {target} from M9 = "
          f"{sol.fobj or 0.0:.4f} ({status})")
    if not sol.fobj or sol.fobj <= MIN_FLUX:
        return 0

    # Every reaction carrying flux that touches ammonium. `nh4` enters the
    # network exactly here, so this is the shortlist of entry points.
    print(f"\n{'reaction':16s} {'flux':>10s}  {'nh4 coeff':>9s}  gene association")
    rows = []
    for rid, flux in sol.values.items():
        if abs(flux) <= MIN_FLUX or rid not in model.reactions:
            continue
        reaction = model.reactions[rid]
        coeff = reaction.stoichiometry.get("M_nh4_c")
        if coeff is None:
            continue
        gpr = reaction.gpr
        genes = str(gpr) if gpr is not None else ""
        rows.append((rid, flux, coeff, genes or "NONE (gap-fill)"))

    for rid, flux, coeff, genes in sorted(rows, key=lambda r: -abs(r[1])):
        print(f"{rid:16s} {flux:10.4f}  {coeff:9.1f}  {genes}")

    if not rows:
        print("  no ammonium-consuming reaction carries flux — "
              "the nitrogen comes in some other way")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
