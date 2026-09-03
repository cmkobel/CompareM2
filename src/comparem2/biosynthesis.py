"""What each genome's metabolic model can build, and what it must acquire.

CarveMe hands the pipeline a network. This reads a high-level phenotype off it:
for each of 32 building blocks — the twenty amino acids, ten vitamins and
cofactors, two quinones — is there a complete, connected route to that compound
from a minimal medium?

**Why not simulate growth on a medium.** Because measured on eleven real draft
models (four *E. faecium*, seven *S. aureus*), **none of them grows on any
defined medium**: M9, M9 anaerobic, M9 glycerol, LB and LB anaerobic all return
exactly zero, and every one of the eleven grows only on the complete medium.
The reason is that growth is a single bit that one unreachable metabolite
destroys: on LB, `116_2` can make 52 of its 53 biomass precursors and fails on
menaquinol-8; `COL` fails on asparagine alone. Per-compound, 53 bits survive
what kills the one. The `media` table below reports that check per genome
rather than hiding it.

**Three verdicts, not two.** A plain producibility scan on M9 cascades: no
folate means no purines means no ATP means everything is blocked, and `116_2`
comes out with 26 of 53 precursors unreachable, which reads as 26 auxotrophies
and is not. So each compound is tried twice —

    de_novo   reachable from M9 alone: salts, glucose, ammonium, phosphate,
              sulfate, oxygen. The genome can build it from scratch.
    upstream  not reachable from M9, but reachable from M9 plus every *other*
              panel compound. The route exists; something else on the panel is
              what is missing.
    none      not reachable even then. No route in this draft model.
    absent    the compound is not in the model at all.

The background for `upstream` is M9 plus the panel, never the complete medium,
and that is load-bearing: on a complete medium `COL` looks able to make
asparagine, because it can take up the Gly-Asn dipeptide and hydrolyse it.
Salvage is not synthesis.

**The probe target is the form the cell needs, not the form on the vitamin
bottle**, and getting that wrong is the failure mode here. Three targets were
tried and rejected against `iML1515`, where the right answer is known:

- `fol` — folate is not an intermediate of de novo synthesis, which runs
  dihydropteroate → dihydrofolate → THF. Probing it called *S. aureus* unable
  to make the compound sulfonamides work by blocking. Target `thf`.
- `thm` — free thiamine is a salvage substrate; de novo synthesis ends at
  thiamine phosphate → ThDP. Probing it called *E. coli* thiamine-auxotrophic.
  Target `thmpp`.
- `lipoate` — free lipoate is salvage too; the de novo product is protein-bound
  lipoyl. Probing it returned "cannot make" for all eleven drafts *and* would
  for any organism. Dropped.

For the same reason no two panel members come from one nutrient family: if both
`nac` and `nad` were on the panel, each would rescue the other and the pair
would report a kinase rather than a pathway.

**Validated against a curated model.** `iML1515` — *E. coli* K-12, manually
curated — returns 31 of 32 de novo. The single exception is adenosylcobalamin,
which *E. coli* genuinely cannot synthesise de novo.

On the drafts the same probe recovers described requirements. All four
*E. faecium* have no route to leucine, methionine, threonine, tryptophan,
valine, riboflavin, pantothenate, NAD and biotin, and three of the four to
arginine and histidine as well; all seven *S. aureus* to thiamine diphosphate,
NAD and biotin. Menaquinone-8 is de novo in four *S. aureus* and unreachable in
every *E. faecium*, and ubiquinone-8 is unreachable or unrepresented in all
eleven — both right for Firmicutes, which use menaquinone and not ubiquinone.

It gets things wrong too, and visibly: all seven *S. aureus* come out with no
route to asparagine, which is not a described requirement of that organism.
Read a verdict as a statement about the draft model.

**It runs in the tool's environment, under a bare `python`**, like
`carve_scip.py` and for the same reason: it imports `reframed`, which comes
with CarveMe and not with CompareM2. Like that module it imports nothing from
its own package, because under `--use-conda` its package is not there at all —
and `reframed` itself is imported inside the functions that need it, so
`report.py` can read `PANEL` from here without needing the solver installed.
A test enforces every part of that.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import NamedTuple

# The paper's phenotype-array protocol: "a maximum uptake rate of 10 mmol/gDW/h
# for every compound" (Machado et al. 2018, Nucleic Acids Research 46:7542).
MAX_UPTAKE = 10.0

# What counts as reachable. The demand reaction's optimum is a feasibility
# witness, not a rate — it saturates at the flux bound — so only its sign is
# read, and this is the tolerance below which a solver's zero is zero.
MIN_FLUX = 1e-6

# Prefixed so it cannot collide with a demand reaction the model already has:
# BiGG ships `DM_4crsol_c` in the E. coli reconstructions.
DEMAND_PREFIX = "R_CM2_DM_"

# M9 and LB as CarveMe defines them, from `carveme/data/input/media_db.tsv`
# (CarveMe, Apache-2.0, © 2017 Daniel Machado). Embedded rather than read out of
# the installed package: these are the media the report names, and a path inside
# another project's package directory is not one to depend on. The anaerobic
# variants are derived by dropping oxygen, which is how media_db defines them.
M9 = ("ca2", "cl", "cobalt2", "cu2", "fe2", "fe3", "glc__D", "h", "h2o", "k",
      "mg2", "mn2", "mobd", "na1", "nh4", "ni2", "o2", "pi", "so4", "zn2")

LB = ("adn", "ala__L", "amp", "arg__L", "aso3", "asp__L", "ca2", "cbl1", "cd2",
      "cl", "cmp", "cobalt2", "cro4", "cu2", "cys__L", "dad_2", "dcyt", "fe2",
      "fe3", "fol", "glc__D", "glu__L", "gly", "gmp", "gsn", "h", "h2o", "h2s",
      "hg2", "his__L", "hxan", "ile__L", "ins", "k", "leu__L", "lipoate",
      "lys__L", "met__L", "mg2", "mn2", "mobd", "na1", "nac", "nh4", "ni2",
      "o2", "phe__L", "pheme", "pi", "pnto__R", "pro__L", "pydx", "ribflv",
      "ser__L", "so4", "thm", "thr__L", "thymd", "trp__L", "tyr__L", "ump",
      "ura", "uri", "val__L", "zn2")

AEROBE = "o2"


class Compound(NamedTuple):
    """One panel entry: the BiGG metabolite probed, and how to name it."""

    bigg: str  # probed as M_<bigg>_c
    name: str
    group: str


AMINO_ACID = "Amino acid"
COFACTOR = "Vitamin or cofactor"
QUINONE = "Quinone"

PANEL = (
    Compound("ala__L", "L-Alanine", AMINO_ACID),
    Compound("arg__L", "L-Arginine", AMINO_ACID),
    Compound("asn__L", "L-Asparagine", AMINO_ACID),
    Compound("asp__L", "L-Aspartate", AMINO_ACID),
    Compound("cys__L", "L-Cysteine", AMINO_ACID),
    Compound("gln__L", "L-Glutamine", AMINO_ACID),
    Compound("glu__L", "L-Glutamate", AMINO_ACID),
    Compound("gly", "Glycine", AMINO_ACID),
    Compound("his__L", "L-Histidine", AMINO_ACID),
    Compound("ile__L", "L-Isoleucine", AMINO_ACID),
    Compound("leu__L", "L-Leucine", AMINO_ACID),
    Compound("lys__L", "L-Lysine", AMINO_ACID),
    Compound("met__L", "L-Methionine", AMINO_ACID),
    Compound("phe__L", "L-Phenylalanine", AMINO_ACID),
    Compound("pro__L", "L-Proline", AMINO_ACID),
    Compound("ser__L", "L-Serine", AMINO_ACID),
    Compound("thr__L", "L-Threonine", AMINO_ACID),
    Compound("trp__L", "L-Tryptophan", AMINO_ACID),
    Compound("tyr__L", "L-Tyrosine", AMINO_ACID),
    Compound("val__L", "L-Valine", AMINO_ACID),
    # The active form in each case — see the module docstring on why the
    # vitamin itself is the wrong target for thiamine, folate and niacin.
    Compound("thmpp", "Thiamine diphosphate (B1)", COFACTOR),
    Compound("ribflv", "Riboflavin (B2)", COFACTOR),
    Compound("nad", "NAD (B3)", COFACTOR),
    Compound("pnto__R", "Pantothenate (B5)", COFACTOR),
    Compound("pydx5p", "Pyridoxal 5'-phosphate (B6)", COFACTOR),
    Compound("btn", "Biotin (B7)", COFACTOR),
    Compound("thf", "Tetrahydrofolate (B9)", COFACTOR),
    Compound("adocbl", "Adenosylcobalamin (B12)", COFACTOR),
    Compound("pheme", "Protoheme", COFACTOR),
    Compound("sheme", "Siroheme", COFACTOR),
    Compound("q8", "Ubiquinone-8", QUINONE),
    Compound("mqn8", "Menaquinone-8", QUINONE),
)

DE_NOVO, UPSTREAM, NO_ROUTE, ABSENT = "de_novo", "upstream", "none", "absent"

PANEL_HEADER = ("compound", "name", "group", "verdict")
MEDIA_HEADER = ("medium", "compounds", "present", "growth")


def _status(solution) -> str:
    """reframed reports status as an enum; take its name either way."""
    return str(getattr(solution.status, "value", solution.status))


def _demand(bigg: str) -> str:
    return f"{DEMAND_PREFIX}{bigg}"


def _exchanges(model) -> dict[str, str]:
    """Compound id to exchange reaction id, for the extracellular exchanges."""
    out = {}
    for rid in model.reactions:
        if rid.startswith("R_EX_") and rid.endswith("_e"):
            out[rid[len("R_EX_"):-len("_e")]] = rid
    return out


def add_demands(model) -> list[Compound]:
    """One drain per panel compound the model carries. Returns those compounds.

    Must run before the solver is built, and `add_reaction_from_str` would
    silently invent the metabolite if it were missing, so membership is checked
    rather than assumed.
    """
    present = []
    for compound in PANEL:
        if f"M_{compound.bigg}_c" not in model.metabolites:
            continue
        rid = _demand(compound.bigg)
        if rid in model.reactions:
            raise SystemExit(f"biosynthesis: {rid} already exists in the model")
        model.add_reaction_from_str(f"{rid}: M_{compound.bigg}_c --> ")
        present.append(compound)
    return present


def medium_constraints(by_compound: dict[str, str], upper: dict[str, float | None],
                       drain_ub: dict[str, float | None], compounds,
                       max_uptake: float = MAX_UPTAKE,
                       open_drain: str | None = None) -> dict:
    """Flux bounds making exactly `compounds` available, and nothing else.

    A pure function of four mappings, so what the solver is asked can be checked
    without a solver — which matters, because a medium that quietly leaves one
    exchange open produces a plausible number and the wrong answer.

    - `by_compound`: compound id to its exchange reaction id
    - `upper`: exchange reaction id to that reaction's own upper bound
    - `drain_ub`: drain reaction id to its own upper bound

    Upper bounds are left as the model has them: this decides what is available
    for uptake, and secretion is not ours to re-decide.

    **Every drain is pinned shut and at most one reopened.** A drain left open
    is a free sink, and a free sink can relieve a steady-state constraint
    elsewhere in the network — which would let one compound's probe change
    another's answer, and let the media table report growth the model cannot
    actually achieve.
    """
    out = {rid: (0.0, upper[rid]) for rid in by_compound.values()}
    out.update({rid: (0.0, 0.0) for rid in drain_ub})
    for cid in compounds:
        rid = by_compound.get(cid)
        if rid is not None:
            out[rid] = (-max_uptake, upper[rid])
    if open_drain is not None:
        out[open_drain] = (0.0, drain_ub[open_drain])
    return out


class _Probe:
    """A model with drains attached, and the solver built over it."""

    def __init__(self, model, max_uptake: float = MAX_UPTAKE):
        from reframed.solvers import solver_instance

        self.model = model
        self.max_uptake = max_uptake
        self.present = add_demands(model)
        self.exchanges = _exchanges(model)
        # After the drains, because the solver is built over the model as it
        # stands and would not know about a reaction added later.
        self.solver = solver_instance(model)
        self._upper = {rid: model.reactions[rid].ub
                       for rid in self.exchanges.values()}
        # Each drain's own upper bound, to restore when it is the one being
        # maximised. Read off the model rather than assumed: it is
        # `add_reaction_from_str` that decides it, and an unbounded reaction is
        # `None` here, not `inf`.
        self._drain_ub = {_demand(c.bigg): model.reactions[_demand(c.bigg)].ub
                          for c in self.present}

    def medium(self, compounds, open_drain: str | None = None) -> dict:
        return medium_constraints(self.exchanges, self._upper, self._drain_ub,
                                  compounds, self.max_uptake, open_drain)

    def maximum(self, reaction: str, constraints: dict) -> float:
        from reframed import FBA

        solution = FBA(self.model, objective={reaction: 1},
                       constraints=constraints, solver=self.solver)
        if _status(solution) != "Optimal":
            return 0.0
        return solution.fobj or 0.0

    def growth(self, constraints: dict) -> float:
        """The model's own objective — biomass — under `constraints`."""
        from reframed import FBA

        solution = FBA(self.model, constraints=constraints, solver=self.solver)
        if _status(solution) != "Optimal":
            return 0.0
        return max(0.0, solution.fobj or 0.0)


def verdicts(probe: _Probe, min_flux: float = MIN_FLUX) -> list[tuple[str, ...]]:
    """One row per panel compound, in panel order."""
    present = {c.bigg for c in probe.present}
    others = [c.bigg for c in probe.present]
    rows = []
    for compound in PANEL:
        if compound.bigg not in present:
            rows.append((compound.bigg, compound.name, compound.group, ABSENT))
            continue
        drain = _demand(compound.bigg)
        if probe.maximum(drain, probe.medium(M9, drain)) > min_flux:
            verdict = DE_NOVO
        else:
            background = list(M9) + [c for c in others if c != compound.bigg]
            verdict = (UPSTREAM
                       if probe.maximum(drain, probe.medium(background, drain)) > min_flux
                       else NO_ROUTE)
        rows.append((compound.bigg, compound.name, compound.group, verdict))
    return rows


def media(probe: _Probe) -> list[tuple[str, ...]]:
    """Growth on each reference medium, with how much of it the model can take.

    `present` is the diagnostic: a medium whose compounds the model has no
    exchange for is not the medium it was asked for. The eleven drafts measured
    carry exchanges for 48 to 51 of LB's 65, against 62 for the curated
    `iML1515`, and the missing ones are the vitamins and nucleosides — which is
    why a Gram-positive draft returns zero on a rich medium.
    """
    anaerobic = tuple(c for c in M9 if c != AEROBE)
    lb_anaerobic = tuple(c for c in LB if c != AEROBE)
    rows = []
    for name, compounds in (("M9", M9), ("M9[-O2]", anaerobic),
                            ("LB", LB), ("LB[-O2]", lb_anaerobic)):
        count = sum(1 for c in compounds if c in probe.exchanges)
        growth = probe.growth(probe.medium(compounds))
        rows.append((name, str(len(compounds)), str(count), f"{growth:.4f}"))
    every = list(probe.exchanges)
    rows.append(("complete", str(len(every)), str(len(every)),
                 f"{probe.growth(probe.medium(every)):.4f}"))
    return rows


def write_tsv(path: Path, header, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="biosynthesis",
        description="Biosynthetic capability from a CarveMe model. "
                    "See module docstring.")
    p.add_argument("--model", type=Path, required=True, help="SBML model to read")
    p.add_argument("--output", type=Path, required=True,
                   help="TSV of per-compound verdicts to write")
    p.add_argument("--media", type=Path, required=True,
                   help="TSV of growth on the reference media to write")
    p.add_argument("--max-uptake", type=float, default=MAX_UPTAKE,
                   help=f"mmol/gDW/h per available compound (default {MAX_UPTAKE})")
    p.add_argument("--min-flux", type=float, default=MIN_FLUX,
                   help=f"below this a demand flux is zero (default {MIN_FLUX})")
    args = p.parse_args(argv)

    from reframed import load_cbmodel

    # flavor='bigg' is what gives the R_/M_ prefixes every id here assumes, and
    # it is what CarveMe writes.
    model = load_cbmodel(str(args.model), flavor="bigg")
    probe = _Probe(model, max_uptake=args.max_uptake)

    # The media table first, so its numbers are read off a model whose drains
    # are shut — which they are in every constraint set, but the ordering says
    # so without the reader having to check.
    write_tsv(args.media, MEDIA_HEADER, media(probe))
    rows = verdicts(probe, min_flux=args.min_flux)
    write_tsv(args.output, PANEL_HEADER, rows)

    counts: dict[str, int] = {}
    for _, _, _, verdict in rows:
        counts[verdict] = counts.get(verdict, 0) + 1
    print(f"biosynthesis: {args.model.name} — "
          + ", ".join(f"{n} {k}" for k, n in sorted(counts.items())),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
