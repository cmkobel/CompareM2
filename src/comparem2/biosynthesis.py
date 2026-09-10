"""What each genome's metabolic model can build, and what it must acquire.

CarveMe hands the pipeline a network. This reads a high-level phenotype off it:
for each of 30 building blocks — the twenty amino acids, nine vitamins and
cofactors, one quinone — is there a complete, connected route to that compound
from a minimal medium?

**Why not simulate growth on a medium.** Because growth is a single bit that one
unreachable metabolite destroys. A CarveMe draft of *B. subtilis* 168 — an
organism that grows on glucose and ammonium — answers **29 of 30 compounds
*de novo* and still returns exactly zero** on M9, M9 anaerobic, LB and LB
anaerobic. On LB, `116_2` can make 52 of its 53 biomass precursors and fails on
menaquinol-8; `COL` fails on asparagine alone. Per-compound, 29 bits survive
what kills the one.

**That zero is a property of the genome in front of it, not of the method.** A
draft of *E. coli* K-12 carved the same way — CarveMe 1.6.6 through
`carve_scip.py`, no gap-filling — grows 0.70 h⁻¹ on M9 and 5.20 on LB. All
eleven Firmicute drafts measured here (four *E. faecium*, seven *S. aureus*)
grow only on the complete medium, which is why the `media` table below reports
the check per genome rather than generalising from it.

**Three verdicts, not two.** A plain producibility scan on M9 cascades: no
folate means no purines means no ATP means everything is blocked, and `116_2`
comes out with 26 of 53 precursors unreachable, which reads as 26 auxotrophies
and is not. So each compound is tried twice —

    de_novo   reachable from M9 alone: salts, glucose, ammonium, phosphate,
              sulfate, oxygen. The genome can build it from scratch.
    upstream  not reachable from M9, but reachable from M9 plus every *other*
              panel compound the model can take up. The route exists; something
              else on the panel is what is missing.
    none      not reachable even then. No route in this draft model.
    absent    the compound is not in the model at all.

**`upstream` says why, never whether.** If a compound were reachable from M9
plus one the genome can already make, it would be reachable from M9 — so the
transitive closure of the panel from the `de_novo` set is the `de_novo` set,
measured and unchanged on thirteen models. On the minimal medium an `upstream`
compound is a requirement too; the verdict separates a compound blocked by its
own missing pathway from one blocked by something else on the list. In `116_2`
all six and in `E8202` all nine are rescued only by glycine, L-serine,
L-threonine or L-methionine, and those models can make none of the four.

"the model can take up" is the caveat on the background: six to eight of the
panel have no `R_EX_*_e` in any model measured — the intracellular cofactors —
so they
cannot act as rescuers. Injecting all of them straight into the cytoplasm
instead moved no verdict in six models, so it is a limit worth naming and not
one that currently costs anything.

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

**Two more targets were rejected on a different test: they never vary.** `btn`
and `q8` were on the panel until 2026-09-10 and came out `none` or `absent` in
**12 of 12 CarveMe drafts across 8 species**, `de novo` only in the curated
`iML1515`. Four of those organisms — *E. coli*, *B. subtilis*, *S. oneidensis*,
*R. solanacearum* — are biotin prototrophs, and *E. coli* and *S. oneidensis*
both use ubiquinone-8, so those are false negatives and not lineage. Traced on
a fresh *E. coli* draft: the biotin chain dies at pimeloyl-CoA (`pmcoa_c`,
which `iML1515` does not contain at all — it reaches biotin by the ACP route),
and no reaction in the draft synthesises `q8`, only the eight quinol oxidations
that cycle it. **The routes are in the database and carving does not keep
them**, so a column that is always the same answer was costing every genome two
false dependencies. Menaquinone-8 stays and does vary.

**Validated at two levels.** The BiGG universe itself — every reaction capped
at ±1000, exchanges added for M9, 25,348 reactions — reaches all 32 compounds
that were on the panel then, so every target here is a representation the
database can actually reach and no verdict below is a badly chosen metabolite.
And the curated model and the drafts now agree: `iML1515` returns **29 of 30**,
and so do drafts of *E. coli*, *B. subtilis*, *P. aeruginosa*, *S. oneidensis*
and *R. solanacearum*. In every one the single miss is adenosylcobalamin —
correct for *E. coli*, which cannot synthesise it de novo, and an honest "not in
this draft" for the rest. **29 of 30 is the calibration**, and it is the same
number whether the model was curated by hand or carved by this pipeline.

On the drafts the probe recovers described requirements. All four *E. faecium*
have no route to leucine, methionine, threonine, tryptophan, valine, riboflavin,
pantothenate and NAD, and three of the four to arginine and histidine as well;
all seven *S. aureus* to thiamine diphosphate and NAD. Menaquinone-8 is de novo
in four *S. aureus* and unreachable in every *E. faecium*, which is a real
lineage difference.

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

# M9's two elemental sources. The panel's `de_novo` verdict is measured on M9,
# so a model that cannot get carbon or nitrogen into its cytoplasm from that
# medium answers `de_novo` to nothing — one hole, a zero in every row — and the
# whole column then describes the model rather than the organism.
#
# **Probed in the cytoplasm, not as a membership test on the exchange set**,
# because the two differ and the difference is not hypothetical. Of the eight
# models in the 2026-09-08 showcase run, four could not reach ammonium and they
# failed in two different places: three carried no `EX_nh4_e` at all, and the
# fourth carried both `EX_nh4_e` and `NH4tex` with no periplasm-to-cytoplasm
# `NH4tpp`. An exchange-set check clears that fourth model, which is as broken
# as the other three — it was the one that reached a certified optimum with the
# most reactions in the set.
#
# Only the M9 family gets this. LB carries nitrogen in every amino acid and
# carbon in most of them, so no single compound is its source and a missing
# `nh4` there is not a defect.
SOURCES = (("glc__D", "carbon"), ("nh4", "nitrogen"))


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
    # vitamin itself is the wrong target for thiamine, folate and lipoate.
    Compound("thmpp", "Thiamine diphosphate (B1)", COFACTOR),
    Compound("ribflv", "Riboflavin (B2)", COFACTOR),
    Compound("nad", "NAD (B3)", COFACTOR),
    Compound("pnto__R", "Pantothenate (B5)", COFACTOR),
    Compound("pydx5p", "Pyridoxal 5'-phosphate (B6)", COFACTOR),
    Compound("thf", "Tetrahydrofolate (B9)", COFACTOR),
    Compound("adocbl", "Adenosylcobalamin (B12)", COFACTOR),
    Compound("pheme", "Protoheme", COFACTOR),
    Compound("sheme", "Siroheme", COFACTOR),
    Compound("mqn8", "Menaquinone-8", QUINONE),
)

DE_NOVO, UPSTREAM, NO_ROUTE, ABSENT = "de_novo", "upstream", "none", "absent"

PANEL_HEADER = ("compound", "name", "group", "verdict")
# `missing` names the medium's compounds the model has no exchange for, and
# `unreachable` the elements of `SOURCES` it cannot get into the cytoplasm.
# Both exist because the `present` count could not tell "no nickel" from "no
# nitrogen source": on 2026-09-08 two models of one strain both read 17 of 20
# for M9 and differed in which three were missing.
# `precursors` and `blocked` answer the question a zero in `growth` raises and
# nothing else here could: which of the things biomass needs is the model
# unable to make. Filled only on a medium that did not grow — 53 solves each,
# against 70 for the whole of the rest of this module.
MEDIA_HEADER = ("medium", "compounds", "present", "growth", "missing",
                "unreachable", "precursors", "blocked")


# ReFramed's `Status` enum, by value: Optimal, Unknown, Suboptimal, Unbounded,
# Infeasible, "Infeasible or Unbounded". Only the first two matter here and
# both are named, because the mapping from status to verdict is the part that
# is easy to get silently wrong.
OPTIMAL, UNBOUNDED = "Optimal", "Unbounded"


def _status(solution) -> str:
    """reframed reports status as an enum; take its name either way."""
    return str(getattr(solution.status, "value", solution.status))


def flux_from(status: str, value: float | None) -> float:
    """What a solve reports as the demand's flux, read the way a verdict needs.

    **`Unbounded` is not zero, it is the opposite of zero**, and treating every
    non-optimal status as "no flux" turns it into `none`. Not a hypothetical:
    pointed at the BiGG universe, whose 25,348 reactions ship bounded at ±inf,
    this probe reported carbon *and* nitrogen unreachable and 29 of 32 panel
    compounds `none` — every one of them an unbounded LP misread. With the
    bounds capped the same model answers 32 of 32.

    Every other non-optimal status is a zero. `Infeasible` genuinely is one;
    `Unknown` and `Suboptimal` are not knowledge, and calling a compound
    producible on a solve that failed would be worse than under-reporting it.
    A pure function, so the mapping can be checked without a solver.
    """
    if status == UNBOUNDED:
        return float("inf")
    if status != OPTIMAL:
        return 0.0
    return value or 0.0


def growth_cell(status: str, value: float | None) -> str:
    """The media table's growth cell: a rate, or the reason there is not one.

    `0.0000` for an infeasible LP is the wrong answer to a different question.
    A model with a maintenance floor its medium cannot pay did not grow slowly,
    it could not be solved — and the complete-medium row is the report's
    control that the model is feasible at all, which a zero there cannot say.
    `iML1515` carries `R_ATPM` at `lb = 6.86`; CarveMe drafts carry `0.0`,
    which is why no run has hit this yet.
    """
    if status != OPTIMAL:
        return status.lower()
    return f"{max(0.0, value or 0.0):.4f}"


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


def add_source_demands(model) -> list[tuple[str, str, str | None]]:
    """One drain per `SOURCES` compound, as (compound, element, reaction).

    Same contract as `add_demands` — must run before the solver is built. A
    compound the model does not carry gets `None` for its reaction rather than
    being dropped, because "the metabolite is not in the network" is one of the
    ways a source is unreachable and has to reach the report as such.
    """
    out: list[tuple[str, str, str | None]] = []
    for bigg, element in SOURCES:
        if f"M_{bigg}_c" not in model.metabolites:
            out.append((bigg, element, None))
            continue
        rid = _demand(bigg)
        if rid in model.reactions:
            raise SystemExit(f"biosynthesis: {rid} already exists in the model")
        model.add_reaction_from_str(f"{rid}: M_{bigg}_c --> ")
        out.append((bigg, element, rid))
    return out


PRECURSOR_PREFIX = "R_CM2_BM_"


def biomass_precursors(model) -> list[str]:
    """The metabolites the model's own objective reaction consumes.

    Which is the question a zero in the growth column raises and the `present`
    count cannot answer: *which* of the things biomass needs is the model
    unable to make? Taken off the objective rather than from a list, because
    CarveMe's biomass composition is the model's, not ours.

    **A drain asks for net production, which is stricter than the biomass
    reaction needs for its maintenance term.** `atp_c` is consumed and `adp_c`
    produced in the same reaction, so growth needs the ATP cycle rather than
    net synthesis of adenosine — and this check asks for the latter. It has
    not misled yet, because a model that cannot make adenosine also fails on
    `datp_c` and `gtp_c`, which biomass genuinely does incorporate. Read a
    blocked `atp_c` alongside those rather than alone.
    """
    objective = [rid for rid, weight in model.get_objective().items() if weight]
    consumed = []
    for rid in objective:
        for met, coefficient in model.reactions[rid].stoichiometry.items():
            if coefficient < 0 and met not in consumed:
                consumed.append(met)
    return consumed


def add_precursor_demands(model) -> dict[str, str]:
    """One drain per biomass precursor, as metabolite id to reaction id.

    Same contract as `add_demands` — before the solver is built. Keyed by the
    full metabolite id because the objective is free to reach outside the
    cytoplasm, even though every model measured stays inside it.
    """
    out = {}
    for met in biomass_precursors(model):
        rid = f"{PRECURSOR_PREFIX}{met.removeprefix('M_')}"
        if rid in model.reactions:
            raise SystemExit(f"biosynthesis: {rid} already exists in the model")
        model.add_reaction_from_str(f"{rid}: {met} --> ")
        out[met] = rid
    return out


def short_metabolite(met: str) -> str:
    """`M_mql8_c` to `mql8`, for a column a reader checks against BiGG.

    The compartment stays on anything that is not cytoplasmic, since then it is
    the interesting part.
    """
    return met.removeprefix("M_").removesuffix("_c")


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
        # Every status this probe's solves came back with, counted. Reported at
        # the end of a run because a non-optimal solve is read as a zero and
        # would otherwise be invisible in a table of zeros.
        self.statuses: dict[str, int] = {}
        self.present = add_demands(model)
        self.sources = add_source_demands(model)
        self.precursors = add_precursor_demands(model)
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
        # The source and precursor drains are pinned shut with the rest, so
        # probing one cannot leave another open as a free sink. The precursor
        # ones matter most here: there are 53 of them on a CarveMe draft, and
        # 53 open sinks would let the media table report growth the model
        # cannot achieve.
        self._drain_ub.update({rid: model.reactions[rid].ub
                               for _, _, rid in self.sources if rid is not None})
        self._drain_ub.update({rid: model.reactions[rid].ub
                               for rid in self.precursors.values()})

    def medium(self, compounds, open_drain: str | None = None) -> dict:
        return medium_constraints(self.exchanges, self._upper, self._drain_ub,
                                  compounds, self.max_uptake, open_drain)

    def _solve(self, objective, constraints: dict):
        """One FBA, with its status counted. Returns (status, objective value)."""
        from reframed import FBA

        solution = FBA(self.model, objective=objective,
                       constraints=constraints, solver=self.solver)
        status = _status(solution)
        self.statuses[status] = self.statuses.get(status, 0) + 1
        return status, solution.fobj

    def maximum(self, reaction: str, constraints: dict) -> float:
        return flux_from(*self._solve({reaction: 1}, constraints))

    def growth(self, constraints: dict) -> tuple[str, float | None]:
        """The model's own objective — biomass — under `constraints`.

        Returns the status too, because `growth_cell` needs it: an infeasible
        LP is not a growth rate of zero.
        """
        return self._solve(None, constraints)


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


def unreachable_sources(probe: _Probe, compounds,
                        min_flux: float = MIN_FLUX) -> list[str]:
    """Which of `SOURCES` the model cannot reach in the cytoplasm on `compounds`.

    The check that has to pass before any verdict on this medium means
    anything about the organism. See `SOURCES` for why it is a flux probe and
    not a lookup in `probe.exchanges`.
    """
    blocked = []
    for _, element, rid in probe.sources:
        if rid is None or probe.maximum(rid, probe.medium(compounds, rid)) <= min_flux:
            blocked.append(element)
    return blocked


def blocked_precursors(probe: _Probe, compounds,
                       min_flux: float = MIN_FLUX) -> list[str]:
    """Which biomass precursors the model cannot produce on `compounds`."""
    return [met for met, rid in probe.precursors.items()
            if probe.maximum(rid, probe.medium(compounds, rid)) <= min_flux]


def media(probe: _Probe, min_flux: float = MIN_FLUX) -> list[tuple[str, ...]]:
    """Growth on each reference medium, and on a zero, what it was short of.

    `present` was the original diagnostic: a medium whose compounds the model
    has no exchange for is not the medium it was asked for. The eleven drafts
    measured carry exchanges for 48 to 51 of LB's 65, against 62 for the curated
    `iML1515`, and the missing ones are the vitamins and nucleosides.

    **It does not explain a zero on a rich medium**, and used to be offered as
    though it did: a *R. solanacearum* draft and `COL` both carry 51 of LB's 65
    and grow 0.7714 and 0.0000. `precursors` and `blocked` are the check that
    does explain it, and they are one metabolite deep — `116_2` reaches 52 of
    its 53 biomass precursors on LB and fails on menaquinol-8, `COL` on
    asparagine alone. Both figures were hand-computed before this was an
    output and are reproduced by it exactly, which is the reason to trust the
    definition in `biomass_precursors`.

    **Only on a medium that did not grow.** 53 solves a medium is eight times
    the rest of this module put together — 6.5 s against 0.8 s on `116_2` —
    and there is nothing to explain about a medium that worked. A model that
    grows on all five pays nothing.

    **A count was not enough**, so `missing` names them and `unreachable` says
    whether what is gone matters. On 2026-09-08 two models of one strain both
    read 17 of 20 for M9 while differing in which three they lacked — and one of
    those three was ammonium, M9's only nitrogen source, which zeroed every one
    of that model's panel verdicts. 17 of 20 looks fine; `nitrogen` does not.
    """
    anaerobic = tuple(c for c in M9 if c != AEROBE)
    lb_anaerobic = tuple(c for c in LB if c != AEROBE)
    every = tuple(probe.exchanges)
    rows = []
    for name, compounds in (("M9", M9), ("M9[-O2]", anaerobic),
                            ("LB", LB), ("LB[-O2]", lb_anaerobic),
                            ("complete", every)):
        absent = tuple(c for c in compounds if c not in probe.exchanges)
        status, value = probe.growth(probe.medium(compounds))
        # M9 family only — see `SOURCES` on why LB gets no such claim, and the
        # complete medium is every exchange the model has rather than a recipe.
        blocked = (unreachable_sources(probe, compounds, min_flux)
                   if name.startswith("M9") else [])
        if status == OPTIMAL and (value or 0.0) > min_flux:
            reached, short = "", ()
        else:
            missed = blocked_precursors(probe, compounds, min_flux)
            reached = f"{len(probe.precursors) - len(missed)}/{len(probe.precursors)}"
            short = tuple(short_metabolite(m) for m in missed)
        rows.append((name, str(len(compounds)),
                     str(len(compounds) - len(absent)),
                     growth_cell(status, value),
                     " ".join(absent), " ".join(blocked),
                     reached, " ".join(short)))
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
    write_tsv(args.media, MEDIA_HEADER, media(probe, min_flux=args.min_flux))
    rows = verdicts(probe, min_flux=args.min_flux)
    write_tsv(args.output, PANEL_HEADER, rows)

    counts: dict[str, int] = {}
    for _, _, _, verdict in rows:
        counts[verdict] = counts.get(verdict, 0) + 1
    print(f"biosynthesis: {args.model.name} — "
          + ", ".join(f"{n} {k}" for k, n in sorted(counts.items())),
          file=sys.stderr)

    # Said out loud, because a status other than `Optimal` becomes a zero in a
    # table already full of them. Every solve on every draft measured so far
    # has been optimal, so this line is normally the one word.
    other = {k: n for k, n in probe.statuses.items() if k != OPTIMAL}
    if other:
        print(f"biosynthesis: {args.model.name} — "
              + f"{sum(other.values())} of {sum(probe.statuses.values())} "
              + "solves were not optimal: "
              + ", ".join(f"{n} {k}" for k, n in sorted(other.items())),
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
