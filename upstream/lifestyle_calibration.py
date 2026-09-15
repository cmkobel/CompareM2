#!/usr/bin/env python3
"""Does a carved model still know how the organism makes a living?

The biosynthesis panel works because it is protected by mass balance: a demand
on `M_trp__L_c` can only carry flux if the atoms came from the medium, so no
spurious cycle can fake it. **A lifestyle is not like that.** "Can this genome
conserve energy from H2 and CO2" is a question about ATP, and ATP is not an
atom — a thermodynamically infeasible cycle answers yes for free. That is the
first reason this is a measurement and not yet a feature.

The second is `btn` and `q8`: routes that the BiGG universe contains and
carving does not keep, wrong in **12 of 12 drafts across 8 species**. Lifestyle
reactions are exactly the rare ones with thin gene evidence in a universe built
mostly from heterotrophs, so the expected failure is a table that says `-`
everywhere. This script is what would show that before anything is built.

**Status: written 2026-09-11, revised 2026-09-15, still not run** — and largely
overtaken. `lifestyle_ceiling.py` ran on 2026-09-15 and found that the
autotrophy, lithotrophy and methanogenesis reactions are in **none** of the five
universes `carve` can load, so those columns cannot be built by any scoring and
the rows below that would have measured them will read `na` or `-` for a reason
that is not about the genome. What is still worth running here is the acceptor
axis. Read `lifestyle.md` first.

Every number below is a prediction or a literature fact, never a measurement,
and the results table at the bottom of this docstring is empty on purpose.

What the 09-15 revision changed, because a later session will otherwise read the
old shape into it: arm three below is new and is the point of the exercise; two
axes that could not be *asked* now say so instead of saying no; *M. barkeri* is
carved from the archaeal universe; and `--faa` scores the ten benchmark genomes.

## What it measures

Fourteen models: seven curated BiGG reconstructions of organisms with known and
different lifestyles, and CarveMe drafts of the same seven organisms. Three
questions, in order of how much they decide:

1. **Does the model make ATP out of nothing?** Every exchange closed to uptake,
   `ATPM` relaxed to a zero floor, maximise an ATP hydrolysis drain. Nonzero
   means the model has a free-energy cycle and **every energy verdict in that
   row is void** — report it, do not score it. This is MEMOTE's test (Lieven
   et al. 2020, Nature Biotechnology 38:272) and CarveMe drafts are not known
   to pass it.
2. **Does the curated model reproduce what the organism is documented to do?**
   That is the answer key checking itself. A curated model is free to omit a
   pathway its organism has, so a miss here is a fact about the model and caps
   what the draft can be scored against.
3. **Does the draft agree with the curated model of the same organism?** This
   is the question. It is model against model, so it is clean — no literature,
   no judgement, just whether carving preserved the network that answers.
4. **Does a column vary across the genomes a user actually has?** Point
   `--model` at the drafts of the 2026-09-11 review set — four *E. faecium*,
   two *S. aureus*, *E. coli*, *M. genitalium*, already carved on thylakoid —
   and see. A column that is the same answer for every genome is not a column,
   which is what took `btn` and `q8` off the biosynthesis panel, and it is the
   question that decides whether any of this belongs in the report. The first
   three can all pass and this one still fail.

## The two probes, and why there are two

Carbon and energy are orthogonal axes, and one score collapses them: a genome
that oxidises H2 for ATP and takes its carbon from acetate is a real lifestyle
(chemolithoheterotrophy) that a growth test scores as a failure of autotrophy.
So each cell gets both —

    energy  max flux through `atp_c + h2o_c --> adp_c + pi_c + h_c`
    carbon  how many of ten central precursors the model can produce

**Both are read as a difference, never as a level.** Each is measured twice:
once on salts + substrate + acceptor, once with the substrate removed. The
verdict credits the substrate only when removing it takes the answer away.

That differential is doing three separate jobs, which is why it is the whole
design rather than a control bolted on:

- it is the local guard against free-energy cycles, since a cycle pays the same
  with the donor and without it;
- it makes the carbon claim honest when the *acceptor* carries carbon. Fumarate,
  DMSO and TMAO are all carbon compounds, so a plain producibility test on
  `glc + dmso` credits glucose for carbon that may have come from DMSO. With
  glucose removed the precursors are still reachable, so no carbon claim is
  made. The honest statement is "the substrate is not necessary", and that is
  what the cell then says;
- it is what lets autotrophy be asked at all. On `h2 + co2`, dropping H2 leaves
  the carbon source in place and the precursors unreachable, which is exactly
  the claim "this genome fixes CO2 and needs the donor to do it".

Verdicts: `CE` substrate supplies both, `E` energy only, `C` carbon only,
`-` neither, `na` the model has no exchange for the substrate or the acceptor
and the question cannot be put to it. `na` is the `absent` of the panel and
matters for the same reason: a strict anaerobe with no `EX_o2_e` has not failed
an aerobic test, it was never asked one.

**And `?`, which is the same distinction one level in.** `na` is about the
medium; `?` is about the model. A model missing any one of the ten precursors
can never satisfy the carbon probe, and one missing a piece of the ATP drain can
never answer the energy probe — in both cases the probe returns a fixed no that
looks exactly like a measurement. It suffixes rather than replaces, so `E?`
reads "energy yes, carbon unaskable". Until 2026-09-15 both cases came out `-`,
which is the conflation `na` exists to prevent, and it pointed the same way as
the prediction on record — the run would have confirmed itself.

## Choices that are not arbitrary

**Ten precursors, not twelve.** The classic set is twelve; `accoa` and `succoa`
are dropped because a drain on them asks for net CoA synthesis — sulfur,
pantothenate, a whole second pathway — rather than for carbon assimilation.
This is the same trap the panel's `atp_c` note describes. The ten that remain
are C/H/O/P only.

**Sulfate and ferric iron come out of the background.** M9 carries `so4` as the
sulfur source and `fe3` as an iron source, which would put the acceptor of a
sulfate reducer and of an iron reducer into every medium and make those two
columns untestable. Neither is needed here: none of the ten precursors contains
sulfur, and `fe2` remains as the iron source. So the salts are M9 minus
glucose, oxygen, sulfate and ferric iron, and `so4` and `fe3` are supplied only
where they are the acceptor being tested.

**`ATPM` is relaxed to a zero floor.** `iML1515` carries `R_ATPM` at `lb=6.86`;
a maintenance floor the medium cannot pay makes the LP infeasible, which is not
the same answer as "no ATP" and would make the gate untestable besides. CarveMe
drafts carry `0.0` already.

**A curated model has a second medium and it is easy to miss.** Closing every
`R_EX_*_e` to uptake is the whole medium for a CarveMe draft, but BiGG models
also ship demand and sink reactions, and a sink with a negative lower bound is
a compound arriving from nowhere. `iJN1463` has one — `R_SK_pqqA_kt_c` at
`lb = -1`, a peptide's worth of carbon — and one is enough to feed a precursor
probe. They are pinned as sources and left alone as sinks.

**No gap-filling *of the draft*, and that is not an oversight.** Gap-fill a model
to a medium and then ask whether it grows on that medium and the answer is yes
by construction. What that argument does not reach is cost-as-score, which is
arm three below: *k* = 0 and *k* = 12 are different answers and neither is
assumed. The 09-11 dismissal was too broad and the 09-15 entry in `DECISIONS.md`
records the reversal.

## Arm three: interrogate the evidence, not the draft

The draft is the wrong object. Carving is a global parsimony MILP, and a rare
pathway with a handful of moderately-scored reactions loses to that objective
routinely — that is what `btn` and `q8` were, wrong in 12 of 12 drafts across 8
species. Asking a draft whether it does Wood–Ljungdahl asks a question the MILP
already answered on unrelated grounds.

`carveme/<sample>.tsv` is the DIAMOND hit table, on disk for every genome this
pipeline has ever run, and `reaction_scoring` turns it into a per-reaction
evidence score over the whole universe — which is exactly what the MILP consumes
and then discards. So `gap_count` restricts nothing away, asks the lifestyle
question on the universe, and reports **how many reactions in the flux solution
this genome has no gene evidence for**. That count is module completeness
generalised from a linear pathway to a network, and it degrades gracefully where
a verdict on a draft collapses to `-`.

It is an **upper bound**: pFBA minimises total flux, not the number of
evidence-free reactions, so *k* = 0 is exact and means the genome's own evidence
suffices, while *k* = 7 does not mean seven are needed. The exact count is a
MILP per cell and is the next thing to build if the bound proves loose. Energy
axis only — see `gap_count`.

## The seven organisms

Chosen because BiGG has a curated model and UniProt or CarveMe has a proteome,
and because between them they cover six ways of making a living. The expected
sets are properties of the **organism**, from the literature, not of the model:

| organism | curated | expected |
| --- | --- | --- |
| *E. coli* K-12 MG1655 | `iML1515` | aerobic, acetate, nitrate, fumarate, DMSO, TMAO, fermentation |
| *B. subtilis* 168 | `iYO844` | aerobic, acetate, nitrate, fermentation |
| *Synechocystis* sp. PCC 6803 | `iJN678` | photoautotrophy, glucose |
| *M. barkeri* | `iAF692` | H2/CO2, CO, methanol and acetate methanogenesis |
| *C. ljungdahlii* DSM 13528 | `iHN637` | H2/CO2, CO, sugar fermentation |
| *G. metallireducens* GS-15 | `iAF987` | Fe(III) on acetate |
| *P. putida* KT2440 | `iJN1463` | aerobic, acetate |

*M. barkeri* is the one imperfect pair: the curated model is str. Fusaro and
UniProt has no Fusaro proteome, so the draft is carved from str. MS
(`UP000033033`). Methanogenesis is conserved well within the species, but the
draft-versus-curated count for that row is strain-to-strain as well as
draft-to-curated, and should be read accordingly.

**Two of the fifteen columns have no positive in the set**, and they are the
specificity check rather than an omission: BiGG has no curated sulfate reducer
this can load and no curated aerobic methylotroph, so every one of the fourteen
models should answer `-` or `na` to `lac_so4` and `meoh_o2`. Anything that
claims either is a false positive worth chasing — a table read only for its
hits cannot tell a working probe from a generous one.

## Running it

Self-contained, like `carve_scip.py`, `biosynthesis.py` and
`panel_calibration.py` — the tool environment's own python, `biosynthesis`
imported by path, nothing from `comparem2`. **Activate the environment, do not
just call its interpreter**: `carve` shells out to DIAMOND.

    conda activate <the carveme env>
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle --curated-only
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle --organism mbar
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle --evidence

    # question 4, on genomes whose lifestyles are known to differ
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle --organism none \
        --faa <benchmark>/*.faa --evidence

    # question 4 again, on a set a real user would have
    python3 lifestyle_calibration.py --workdir /tmp/lifestyle --organism none \
        --model ~ghrunner/cm2-report-run/*/carveme/*.xml

The last two are different questions and both are needed. The benchmark genomes
ask whether the probe can discriminate *at all*; a review set of four
*E. faecium*, two *S. aureus*, *E. coli* and *M. genitalium* asks whether it
discriminates on the genomes CompareM2's users bring. A flat column on the
second set alone cannot tell a blind probe from a set of similar organisms.

Downloads are cached in the workdir and reused. Models come from
`bigg.ucsd.edu/static/models/`, proteomes from UniProt's REST stream, and
*E. coli* and *B. subtilis* from the proteomes CarveMe already bundles. All
nine identifiers were checked against the live endpoints on 2026-09-11.

Budget: about 350 solves a model — 15 lifestyles by (3 energy + up to 20
carbon) — so fourteen models is minutes, not hours. Carving the five downloaded
proteomes is the slow half and happens once. `--evidence` is the expensive one
and is opt-in for that reason: it loads the whole carving universe per genome
and solves over it, against a draft's ~1,500 reactions. That universe is
**5,532 reactions, not 25,348** — the larger figure is `bigg_universe.xml.gz`,
an input to CarveMe's build that `carve` never loads. Measured 2026-09-15.

## The benchmark set, and a defect it exposed before it was run

`BENCHMARK` is ten complete RefSeq genomes, one per way of making a living,
`--benchmark` fetches and scores them, and every accession was resolved against
the NCBI API rather than written from memory. It exists because `ORGANISMS` has
no positive at all for `lac_so4` or `meoh_o2` — those two columns could only
ever be confirmed negative — and because question four needs a set where a
column *must* vary if the probe works.

**And it found this without a solver running.** *S. oneidensis* MR-1 is the
organism that respires nitrate, fumarate, DMSO, TMAO and Fe(III) — the whole
acceptor axis, which is the only axis expected to survive gate 4 — and it does
not catabolise glucose. Every acceptor column here pairs its acceptor with
glucose, so MR-1 answers `-` down that axis for a reason that is about the donor
and says nothing about the acceptor. **A donor baked into the column is a design
defect**, and the shape of it is: `glc_fum`, `glc_dmso` and `glc_tmao` have
exactly one organism expected to fire them (*E. coli*), and the natural second
positive for all three is the organism the donor excludes. Whether to fix it by
standardising on a more widely catabolised donor, or by crediting an acceptor
when any of a small donor panel works, is a design decision and is not taken
here.

## Results

Not yet run — no Linux host was reachable on 2026-09-15. Fill in from one, and
say the date and the CarveMe and ReFramed versions, as `panel_calibration.py`
does.

    gate failures         —
    could not be gated    —
    curated vs literature —
    draft vs curated      —
    benchmark vs literature —
    columns that vary     —
    evidence-restricted k —
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import math
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import NamedTuple

HERE = Path(__file__).resolve().parent
BIOSYNTHESIS = HERE.parent / "src" / "comparem2" / "biosynthesis.py"
CARVE_SCIP = HERE.parent / "src" / "comparem2" / "carve_scip.py"

BIGG_URL = "https://bigg.ucsd.edu/static/models/{model}.xml.gz"
UNIPROT_URL = ("https://rest.uniprot.org/uniprotkb/stream"
               "?query=proteome:{upid}&format=fasta")

# Verdicts. `na` is the panel's `absent` — the model has no exchange for this
# substrate or acceptor, so the question could not be put to it.
BOTH, ENERGY, CARBON, NEITHER, NOT_TESTABLE = "CE", "E", "C", "-", "na"

# `?` is the *other* way a question goes unasked, and it is a property of the
# model rather than of the cell: a model with no `M_e4p_c` can never satisfy the
# carbon probe, and one missing a piece of the ATP drain can never answer the
# energy probe. It suffixes rather than replaces, so a claim that was actually
# established still reads — `E?` is "energy yes, carbon unaskable".
#
# Without it both cases came out `-`, which turns "could not ask" into "the
# answer is no" — the exact conflation `na` exists to prevent, and in the
# direction of the prediction on record, so the run would have confirmed itself.
UNASKABLE = "?"

# The ten central precursors, C/H/O/P only. `accoa` and `succoa` are off the
# classic twelve on purpose — see the module docstring.
PRECURSORS = ("g6p", "f6p", "r5p", "e4p", "g3p", "3pg", "pep", "pyr", "oaa",
              "akg")

CARBON_PREFIX = "R_CM2_C_"
ATP_DRAIN = "R_CM2_DM_atp"

ATP_DRAIN_EQUATION = (f"{ATP_DRAIN}: M_atp_c + M_h2o_c --> "
                      "M_adp_c + M_pi_c + M_h_c")

# Every metabolite that equation needs, checked rather than assumed. ReFramed's
# `add_reaction_from_str` *creates* a metabolite it does not recognise, so a
# model missing `M_pi_c` would get a drain onto a compound nothing produces and
# report no ATP — a model defect read as a biological answer.
ATP_DRAIN_PARTS = ("atp", "h2o", "adp", "pi", "h")

# Anything the model's own maintenance reaction forces. Relaxed to a zero floor
# before the solver is built, because a floor the medium cannot pay is an
# infeasible LP and not an answer.
MAINTENANCE = ("R_ATPM", "R_NGAM", "R_ATPM_c")

# Demand and sink reactions the model already carries. Blocked as *sources* —
# see `block_sources`.
SOURCE_PREFIXES = ("R_DM_", "R_SK_", "R_sink_")

# A reaction counts as evidence-backed when this genome's DIAMOND hits give it a
# positive normalised score. Zero rather than something tuned: CarveMe splits
# there itself (`default_score=-1.0` for everything unscored), and what a better
# cutoff would be is unmeasured — it is one of the plan's open questions.
EVIDENCE_MIN = 0.0

# What arm three says when there is no number to give.
# How much more ATP the acceptor has to buy before the column calls it
# respiration. **Uncalibrated**: it sits in the gap between `iYO844`'s 1.12 on
# fumarate, which is a false positive, and `iML1515`'s 1.40, which is true — two
# points, so this is a placeholder with a reason, not a measurement. See
# `Probe.respires`.
YIELD_MIN = 1.2

NO_ROUTE = "none"       # not even the whole universe can do it: a database limit
NOT_NEEDED = "free"     # it happens without the substrate, so no claim is made


class Lifestyle(NamedTuple):
    """One way of making a living, as a medium and what to take out of it.

    `substrate` is the electron donor and candidate carbon source; it is what
    the differential removes. `acceptor` is the terminal acceptor, present in
    both media, and empty for a fermentation or a methanogenesis where the
    donor's own carbon is the sink.
    """

    key: str
    name: str
    substrate: tuple[str, ...]
    acceptor: tuple[str, ...]
    # Whether `acceptor` is a *terminal electron acceptor*, which decides how the
    # cell is scored. The field is otherwise overloaded: on `h2_co2`, `co_co2`
    # and `photo_co2` it holds CO2, which is the carbon source rather than
    # somewhere electrons are put. A yield ratio is the right question for the
    # first kind and meaningless for the second — light alone already makes ATP,
    # so photoautotrophy would score a ratio of 1 and read negative.
    respiratory: bool = False


LIFESTYLES = (
    Lifestyle("glc_o2", "Aerobic heterotrophy, glucose", ("glc__D",), ("o2",),
              respiratory=True),
    Lifestyle("ac_o2", "Acetate oxidation, aerobic", ("ac",), ("o2",),
              respiratory=True),
    Lifestyle("glc_no3", "Nitrate respiration", ("glc__D",), ("no3",),
              respiratory=True),
    Lifestyle("glc_fum", "Fumarate respiration", ("glc__D",), ("fum",),
              respiratory=True),
    Lifestyle("glc_dmso", "DMSO respiration", ("glc__D",), ("dmso",),
              respiratory=True),
    Lifestyle("glc_tmao", "TMAO respiration", ("glc__D",), ("tmao",),
              respiratory=True),
    Lifestyle("glc_ferm", "Fermentation, glucose", ("glc__D",), ()),
    Lifestyle("ac_fe3", "Fe(III) respiration, acetate", ("ac",), ("fe3",),
              respiratory=True),
    Lifestyle("lac_so4", "Sulfate respiration, lactate", ("lac__D",), ("so4",),
              respiratory=True),
    Lifestyle("h2_co2", "Hydrogenotrophy", ("h2",), ("co2",)),
    Lifestyle("co_co2", "Carboxydotrophy", ("co",), ("co2",)),
    Lifestyle("photo_co2", "Photoautotrophy", ("photon",), ("co2",)),
    Lifestyle("meoh_o2", "Methylotrophy, aerobic", ("meoh",), ("o2",),
              respiratory=True),
    Lifestyle("meoh_anox", "Methylotrophy, anoxic", ("meoh",), ()),
    Lifestyle("ac_anox", "Acetate, no acceptor", ("ac",), ()),
)

LIFESTYLE_KEYS = tuple(life.key for life in LIFESTYLES)


class Organism(NamedTuple):
    """A curated model, a proteome to carve, and what the organism is known for.

    `expected` is literature about the organism and deliberately not about
    either model. `carveme` names a proteome CarveMe bundles; `uniprot` a
    proteome to fetch. Exactly one of the two is set.
    """

    key: str
    name: str
    bigg: str
    carveme: str | None
    uniprot: str | None
    expected: frozenset[str]
    # Which pre-built CarveMe universe to carve from. `None` is the default
    # bacterial one, which is also what `catalogue.py` gives every genome — so
    # anything set here is a deviation from what the pipeline does today, and
    # says so at the point it deviates.
    universe: str | None = None


ORGANISMS = (
    Organism("ecoli", "E. coli K-12 MG1655", "iML1515",
             "Ecoli_K12_MG1655", None,
             frozenset({"glc_o2", "ac_o2", "glc_no3", "glc_fum", "glc_dmso",
                        "glc_tmao", "glc_ferm"})),
    Organism("bsub", "B. subtilis 168", "iYO844",
             "Bsubtilis_168", None,
             frozenset({"glc_o2", "ac_o2", "glc_no3", "glc_ferm"})),
    Organism("syn", "Synechocystis PCC 6803", "iJN678",
             None, "UP000001425",
             frozenset({"photo_co2", "glc_o2"})),
    # The curated model is str. Fusaro and UniProt has no Fusaro proteome, so
    # the draft is str. MS. See the docstring: this pair is strain-to-strain.
    #
    # **Carved from the archaeal universe**, because `MCR` and `HDR` live in
    # `iAF692` — an archaeal model — and so can only be in `universe_archaea`.
    # From the default bacterial universe methanogenesis is unreachable by
    # construction and this row would measure that rather than the organism.
    # `carve`'s flag is verified (`carve.py:280`, `-u/--universe`) and
    # `universe_archaea.xml.gz` is one of the five CarveMe ships; in the
    # pipeline the same thing is spelled `--set carveme--universe=archaea`.
    Organism("mbar", "M. barkeri (curated Fusaro, draft MS)", "iAF692",
             None, "UP000033033",
             frozenset({"h2_co2", "co_co2", "meoh_anox", "ac_anox"}),
             universe="archaea"),
    Organism("clju", "C. ljungdahlii DSM 13528", "iHN637",
             None, "UP000077020",
             frozenset({"h2_co2", "co_co2", "glc_ferm"})),
    Organism("gmet", "G. metallireducens GS-15", "iAF987",
             None, "UP000007073",
             frozenset({"ac_fe3"})),
    Organism("pput", "P. putida KT2440", "iJN1463",
             None, "UP000000556",
             frozenset({"glc_o2", "ac_o2"})),
)


class Benchmark(NamedTuple):
    """A complete genome with a documented lifestyle and no curated model.

    Separate from `ORGANISMS` because there is no answer key: nothing to compare
    a draft against but the literature. What these are for is question four —
    whether a column *varies* — on a set chosen so that it must. A review set of
    four *E. faecium*, two *S. aureus*, *E. coli* and *M. genitalium* cannot
    distinguish a probe that does not discriminate from a set of organisms that
    are alike, and these ten can.
    """

    key: str
    name: str
    accession: str
    expected: frozenset[str]
    universe: str | None = None


# Ten complete genomes, one per way of making a living. Every accession was
# resolved against the NCBI datasets API on 2026-09-15 and none was written from
# memory — two would have been wrong: *G. metallireducens* GS-15 is
# `GCF_000012925.1`, and RefSeq has no *M. barkeri* Fusaro at all, so the
# reference is str. **MS**, which is the strain this file already carves from
# UniProt (`UP000033033`).
#
# The last three have no curated BiGG model and are here because of that:
# `lac_so4` and `meoh_o2` have no positive anywhere in `ORGANISMS`, so without
# `dvul` and `mext` they are pure specificity columns with nothing to confirm
# they can ever fire.
BENCHMARK = (
    Benchmark("ecoli", "E. coli K-12 MG1655", "GCF_000005845.2",
              frozenset({"glc_o2", "ac_o2", "glc_no3", "glc_fum", "glc_dmso",
                         "glc_tmao", "glc_ferm"})),
    Benchmark("bsub", "B. subtilis 168", "GCF_000009045.1",
              frozenset({"glc_o2", "ac_o2", "glc_no3", "glc_ferm"})),
    Benchmark("pput", "P. putida KT2440", "GCF_000007565.2",
              frozenset({"glc_o2", "ac_o2"})),
    Benchmark("syn", "Synechocystis sp. PCC 6803", "GCF_000009725.1",
              frozenset({"photo_co2", "glc_o2"})),
    Benchmark("clju", "C. ljungdahlii DSM 13528", "GCF_000143685.1",
              frozenset({"h2_co2", "co_co2", "glc_ferm"})),
    Benchmark("gmet", "G. metallireducens GS-15", "GCF_000012925.1",
              frozenset({"ac_fe3"})),
    Benchmark("mbar", "M. barkeri MS", "GCF_000970025.1",
              frozenset({"h2_co2", "co_co2", "meoh_anox", "ac_anox"}),
              universe="archaea"),
    # Sulfate reducer. The one organism that can make `lac_so4` fire.
    Benchmark("dvul", "N. (D.) vulgaris Hildenborough", "GCF_000195755.1",
              frozenset({"lac_so4"})),
    # Aerobic methylotroph, and it also grows on C2/C4 acids. The one organism
    # that can make `meoh_o2` fire.
    Benchmark("mext", "M. extorquens AM1", "GCF_000022685.1",
              frozenset({"meoh_o2", "ac_o2"})),
    # **This one is here to break the panel, and it should.** MR-1 respires
    # nitrate, fumarate, DMSO, TMAO and Fe(III) — the whole acceptor axis, the
    # only axis expected to survive gate 4 — and it **cannot catabolise
    # glucose**. Every acceptor column pairs its acceptor with glucose, so MR-1
    # answers `-` to all of them for a reason that is about the donor and
    # nothing to do with the acceptor. Recorded as `ac_o2` only, which is what
    # the columns as written can legitimately ask it. See the module docstring:
    # a donor baked into the column is a design defect the benchmark surfaced
    # before any compute was spent on it.
    Benchmark("sone", "S. oneidensis MR-1", "GCF_000146165.2",
              frozenset({"ac_o2"})),
)

NCBI_URL = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
            "{accession}/download?include_annotation_type=PROT_FASTA")


def benchmark_proteome(entry: Benchmark, workdir: Path) -> Path:
    """RefSeq's proteome for a benchmark genome, downloaded once.

    RefSeq rather than UniProt so all ten come from one annotation pipeline —
    and for the seven that are also in `ORGANISMS`, scoring both is a free
    control on whether the proteome source changes the verdict.
    """
    faa = workdir / f"{entry.key}.faa"
    if faa.exists() and faa.stat().st_size:
        return faa
    archive = fetch(NCBI_URL.format(accession=entry.accession),
                    workdir / f"{entry.key}.zip")
    with zipfile.ZipFile(archive) as bundle:
        name = next(n for n in bundle.namelist() if n.endswith("protein.faa"))
        faa.write_bytes(bundle.read(name))
    archive.unlink()
    return faa


def load_biosynthesis():
    """By path, because under `--use-conda` `comparem2` is not installed here."""
    spec = importlib.util.spec_from_file_location("biosynthesis", BIOSYNTHESIS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def salts(bio) -> tuple[str, ...]:
    """M9 without the things this asks about.

    Glucose and oxygen are substrate and acceptor here; `so4` and `fe3` are the
    acceptors of a sulfate reducer and an iron reducer and cannot also be
    background. `fe2` stays as the iron source, and none of `PRECURSORS`
    contains sulfur.
    """
    out = tuple(c for c in bio.M9 if c not in ("glc__D", "o2", "so4", "fe3"))
    assert "fe2" in out and "nh4" in out, out
    return out


def block_sources(model) -> list[str]:
    """Stop the model's own sinks from running backwards as free supply.

    `medium_constraints` closes every `R_EX_*_e` to uptake, which is the whole
    medium for a CarveMe draft. A curated model also carries demand and sink
    reactions, and a sink with a negative lower bound is a compound arriving
    from nowhere: `iJN1463` ships `R_SK_pqqA_kt_c` at `lb = -1`, a peptide's
    worth of carbon a precursor probe would happily assimilate. One of 51 such
    reactions across the seven curated models, and one is enough.

    Their upper bounds are left alone. Removing a dead-end byproduct is what
    they are for and models are infeasible without them; supplying one is not.
    """
    blocked = []
    for rid in model.reactions:
        if rid.startswith(SOURCE_PREFIXES) and model.reactions[rid].lb < 0:
            model.reactions[rid].lb = 0.0
            blocked.append(rid)
    return blocked


def fetch(url: str, path: Path, unzip: bool = False) -> Path:
    """Download once. An empty body is an error, not an empty file.

    UniProt answers 200 with zero bytes for a proteome that has been superseded
    — `UP000001656` does, which is why *C. ljungdahlii* is carved from
    `UP000077020`. Writing that as a file would carve an empty proteome and
    report a model with no lifestyle at all.
    """
    if path.exists() and path.stat().st_size:
        return path
    print(f"  fetching {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as response:
        body = response.read()
    if unzip:
        body = gzip.decompress(body)
    if not body:
        raise SystemExit(f"lifestyle: {url} returned an empty body")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(body)
    return path


def bundled_fasta() -> Path:
    """Where CarveMe keeps the benchmark proteomes in the installed package."""
    import carveme

    return Path(carveme.__file__).parent / "data" / "benchmark" / "fasta"


def proteome(organism: Organism, workdir: Path) -> Path:
    """Where this organism's protein FASTA comes from — bundled, or UniProt."""
    if organism.carveme is not None:
        faa = bundled_fasta() / f"{organism.carveme}.faa"
        if not faa.exists():
            raise SystemExit(f"lifestyle: no bundled proteome {faa}")
        return faa
    return fetch(UNIPROT_URL.format(upid=organism.uniprot),
                 workdir / f"{organism.key}.faa")


class Draft(NamedTuple):
    """What carving leaves behind: the model, and the evidence it was built on."""

    model: Path
    hits: Path


def draft(faa: Path, output: Path, universe: str | None = None) -> Draft:
    """Carve `faa`, through `carve_scip.py` so it is the pipeline's own path.

    Both outputs are returned because the second one is arm three's entire
    input. `carve` derives its DIAMOND output from the *input* path — `carve.py`
    does `blast_output = os.path.splitext(inputfile)[0] + '.tsv'` — and
    `carve_scip.link_input` links the proteome next to the model, so the hits
    land beside it under the proteome's name. That is the same mechanism that
    once overwrote Bakta's feature table, read here as a feature.
    """
    hits = output.parent / f"{faa.stem}.tsv"
    if output.exists() and hits.exists():
        print(f"  {output.name} exists, reusing", file=sys.stderr)
        return Draft(output, hits)
    command = [sys.executable, str(CARVE_SCIP),
               "--faa", str(faa), "--output", str(output)]
    if universe:
        command += ["--universe", universe]
    subprocess.run(command, check=True)
    # `check=True` is not enough. Measured 2026-09-15: run with only the
    # environment's interpreter and `carve` prints "Unable to run diamond" and
    # **exits 0**, so the failure arrived here as a missing file and surfaced
    # four frames later as `IOError: Model file was not found` from the SBML
    # reader. That is the documented trap — activate the environment, do not
    # just call its python — and it deserves to say so where it happens.
    missing = [p for p in (output, hits) if not p.exists()]
    if missing:
        raise SystemExit(
            f"lifestyle: carve exited 0 but wrote no {' and no '.join(p.name for p in missing)}. "
            "The usual cause is DIAMOND not being on PATH: `carve` shells out to "
            "it, so the environment has to be *activated*, not just its "
            "interpreter called.")
    return Draft(output, hits)


class Probe:
    """A model with the drains attached and the solver built over it.

    Deliberately not `biosynthesis._Probe`: that one adds the 30 panel demands,
    the source demands and one drain per biomass precursor — 80-odd reactions
    this does not use. It reuses that module's `medium_constraints`,
    `flux_from` and `MIN_FLUX`, which are the parts that were hard to get right.
    """

    def __init__(self, bio, model, max_uptake: float | None = None):
        from reframed.solvers import solver_instance

        self.bio = bio
        self.model = model
        self.max_uptake = max_uptake or bio.MAX_UPTAKE
        self.statuses: dict[str, int] = {}

        # Before the solver is built, or the solver does not see any of it.
        self.maintenance = [rid for rid in MAINTENANCE if rid in model.reactions]
        for rid in self.maintenance:
            model.reactions[rid].lb = 0.0
        self.blocked = block_sources(model)

        # What the model lacks decides which probe can be *asked* of it, and
        # that has to be recorded rather than folded silently into a negative.
        self.drains: dict[str, str] = {}
        self.missing_atp = tuple(m for m in ATP_DRAIN_PARTS
                                 if f"M_{m}_c" not in model.metabolites)
        if not self.missing_atp:
            model.add_reaction_from_str(ATP_DRAIN_EQUATION)
            self.drains["atp"] = ATP_DRAIN
        self.missing_precursors = tuple(m for m in PRECURSORS
                                        if f"M_{m}_c" not in model.metabolites)
        for met in PRECURSORS:
            if met in self.missing_precursors:
                continue
            rid = f"{CARBON_PREFIX}{met}"
            model.add_reaction_from_str(f"{rid}: M_{met}_c --> ")
            self.drains[met] = rid

        self.exchanges = bio._exchanges(model)
        self.solver = solver_instance(model)
        self._upper = {rid: model.reactions[rid].ub
                       for rid in self.exchanges.values()}
        self._drain_ub = {rid: model.reactions[rid].ub
                          for rid in self.drains.values()}

    @property
    def energy_askable(self) -> bool:
        """False when the model has no ATP drain to maximise."""
        return "atp" in self.drains

    @property
    def carbon_askable(self) -> bool:
        """False when any one of the ten precursors is absent from the model.

        All ten, not most. The carbon probe reads "every precursor reachable",
        so a model short one of them can never satisfy it and every carbon cell
        in its row would be a fixed no rather than a measurement.
        """
        return not self.missing_precursors

    def light(self) -> tuple[str, ...]:
        """Whatever this model calls a photon exchange, by any of its names.

        `iJN678` spells it `EX_photon_e`; other reconstructions carry
        wavelength-resolved ones. Matched rather than assumed, so a model that
        splits light into bands is still asked the question.
        """
        return tuple(c for c in self.exchanges if "photon" in c)

    def has(self, compounds) -> bool:
        return all(c in self.exchanges for c in compounds)

    def medium(self, compounds, open_drain: str | None = None) -> dict:
        return self.bio.medium_constraints(self.exchanges, self._upper,
                                           self._drain_ub, compounds,
                                           self.max_uptake, open_drain)

    def maximum(self, drain: str, compounds) -> float:
        from reframed import FBA

        solution = FBA(self.model, objective={drain: 1},
                       constraints=self.medium(compounds, drain),
                       solver=self.solver)
        status = self.bio._status(solution)
        self.statuses[status] = self.statuses.get(status, 0) + 1
        return self.bio.flux_from(status, solution.fobj)

    def makes_atp(self, compounds) -> bool:
        drain = self.drains.get("atp")
        if drain is None:
            return False
        return self.maximum(drain, compounds) > self.bio.MIN_FLUX

    def respires(self, life: "Lifestyle", background, substrate) -> float:
        """ATP yield with the acceptor over ATP yield without it.

        **Respiration is a yield, not a feasibility**, and that is why neither
        binary differential can score the acceptor axis. Measured on the curated
        models, 2026-09-15:

        - Removing the *donor* asks whether the donor is necessary. Every
          fermenter then reads positive on every acceptor column — `iML1515`
          reports `lac_so4` because lactate is necessary, though sulfate does
          nothing, and E. coli comes out a sulfate reducer.
        - Removing the *acceptor* asks whether the acceptor is necessary. Every
          fermenter now reads negative on every acceptor column, because it can
          always make some ATP without one: the same `iML1515` stops being an
          aerobe.

        E. coli ferments, and it is the textbook organism for the whole acceptor
        axis, so both binary tests fail on the case the axis exists for. The
        ratio separates them cleanly — glucose plus each acceptor against
        glucose alone, uptake capped at `MAX_UPTAKE`:

            glc_o2  2.80    glc_no3 3.00    glc_dmso 1.50
            glc_fum 1.40    glc_tmao 1.50   lac_so4  **1.00**

        Seven of seven right, including the one true negative at exactly 1.00.

        Two caveats that decide how far this can be read. It is a **within-model**
        ratio: two genomes with different ETC stoichiometry are not comparable on
        it, which is what the plan guessed and this confirms. And the threshold
        is **uncalibrated** — `iYO844` returns 1.12 on fumarate, which is a false
        positive, against `iML1515`'s 1.40 which is true, so the default sits in
        a gap 0.28 wide established on two points. Calibrating it needs many more
        models than seven.

        Infinity when the donor alone yields nothing: the acceptor is then not
        merely improving the yield, it is the reason there is one.
        """
        drain = self.drains["atp"]
        alone = self.maximum(drain, tuple(background) + tuple(substrate))
        if alone <= self.bio.MIN_FLUX:
            return math.inf
        full = self.maximum(drain,
                            tuple(background) + tuple(substrate) + life.acceptor)
        return full / alone

    def precursors(self, compounds) -> int:
        """How many of `PRECURSORS` the model can produce net on this medium."""
        return sum(1 for met, rid in self.drains.items()
                   if met != "atp"
                   and self.maximum(rid, compounds) > self.bio.MIN_FLUX)


def free_energy(probe: Probe) -> bool | None:
    """MEMOTE's test: ATP with every exchange closed to uptake.

    True means the model generates free energy and every energy verdict in its
    row is void. **None means the gate could not be run at all** — no ATP drain
    — which is not the same as passing it, and previously read as `ok`.

    Secretion is left open, which is how MEMOTE poses it: the question is
    whether the network can pay for ATP with nothing coming in, not whether it
    can hold its breath.
    """
    if not probe.energy_askable:
        return None
    return probe.makes_atp(())


def verdict(probe: Probe, life: Lifestyle, background,
            yield_min: float = YIELD_MIN) -> str:
    """One cell: does this substrate supply carbon, energy, both or neither.

    Every claim is the difference between the full medium and the same medium
    with the substrate removed. See the module docstring for the three separate
    things that difference is doing.
    """
    substrate = probe.light() if life.key == "photo_co2" else life.substrate
    if not substrate or not probe.has(substrate) or not probe.has(life.acceptor):
        return NOT_TESTABLE

    full = tuple(background) + tuple(substrate) + life.acceptor
    without = tuple(background) + life.acceptor

    energy = carbon = False
    if probe.energy_askable:
        if life.respiratory:
            # Yield alone, and deliberately *not* also the donor differential.
            # Requiring the donor to be necessary vetoes the textbook case:
            # fumarate is a carbon source as well as an acceptor, so `iML1515`
            # makes ATP from fumarate with no glucose, the donor test fails, and
            # E. coli stops being a fumarate respirer at a yield of 1.40. A
            # free-energy cycle contributes to both solves and so damps the
            # ratio rather than inflating it, and the gate catches the models
            # that have one.
            energy = probe.respires(life, background, substrate) >= yield_min
        else:
            energy = probe.makes_atp(full) and not probe.makes_atp(without)

    claimed = (CARBON if carbon else "") + (ENERGY if energy else "")
    if probe.energy_askable and probe.carbon_askable:
        return claimed or NEITHER
    # One axis was never asked, so silence on the other is not a negative.
    return claimed + UNASKABLE


def score(probe: Probe, background,
          yield_min: float = YIELD_MIN) -> dict[str, str]:
    return {life.key: verdict(probe, life, background, yield_min)
            for life in LIFESTYLES}


def positives(row: dict[str, str]) -> set[str]:
    """Which lifestyles this row claims. An energy claim is the claim.

    Read as a substring rather than by set membership, so a partially answered
    cell still counts: `E?` is an energy claim with the carbon axis unasked, and
    it is as much a claim as `E`.
    """
    return {key for key, value in row.items() if ENERGY in value}


def evidence(hits: Path, universe_name: str | None):
    """CarveMe's universe, scored by one genome's own DIAMOND hits.

    This is the object the 2026-09-15 revision says to interrogate *instead of*
    the draft. `reaction_scoring` is what the carving MILP consumes and then
    throws away, so the scores are upstream of every choice the MILP made on
    global-parsimony grounds — which is what loses a rare pathway. Follows the
    recipe in `upstream/README.md`, the in-repo precedent for driving it
    standalone; that it works off a saved `hits.tsv` is still unverified.

    Returns the universe, the per-reaction score, and the set of reactions BiGG
    knows *any* GPR for — because "this genome has no evidence" and "no gene
    could ever be evidence" are different failures. A transporter in the flux
    solution should not be counted against the genome like a missing enzyme.
    """
    import pandas as pd
    from carveme import config, project_dir
    from carveme.reconstruction.diamond import load_diamond_results
    from carveme.reconstruction.scoring import reaction_scoring
    from reframed import load_cbmodel
    from reframed.core.transformation import apply_bounds

    generated = project_dir + config.get("generated", "folder")
    path = (f"{generated}universe_{universe_name}.xml.gz" if universe_name
            else project_dir + config.get("generated", "default_universe"))
    universe = load_cbmodel(path, flavor="bigg")
    universe.id = "evidence"
    # The universe ships every reaction at ±inf. Leaving that is what made
    # `universe_ceiling.py` read `Unbounded` as zero and report 29 of 32 `none`.
    apply_bounds(universe)
    gprs = pd.read_csv(project_dir + config.get("generated", "bigg_gprs"))
    gprs = gprs[gprs.reaction.isin(list(universe.reactions))]
    scores, _ = reaction_scoring(load_diamond_results(str(hits)), gprs)
    return (universe,
            dict(scores[["reaction", "normalized_score"]].values),
            set(gprs.reaction))


def gap_count(probe: Probe, life: Lifestyle, background, scores, associable) -> str:
    """How many evidence-free reactions the energy route still needs: *k*.

    The differential is kept exactly as `verdict` keeps it — a route that works
    without the substrate makes no claim about the substrate — and here it is
    not optional: the universe certainly contains free-energy cycles, being
    every reaction BiGG has, and the difference is the only thing that cancels
    them.

    **The count is an upper bound, not the minimum.** pFBA minimises total flux,
    not the number of evidence-free reactions; the true minimum is a MILP per
    cell, which is the gap-fill-distance formulation the 09-15 entry
    rehabilitated and the next thing to build if this bound proves loose.
    Saying "upper bound" is load-bearing — *k* = 0 is exact and means the
    genome's own evidence suffices, but *k* = 7 does not mean seven are needed.

    Reported `k+t`: *k* reactions BiGG knows a gene for and this genome does not
    have, *t* reactions no gene could ever explain — transport, spontaneous.
    The plan's open question is what weight those deserve; keeping them apart is
    the least this can do while that is unmeasured. Exchanges are excluded
    outright: they are the medium, which was chosen, not a gap.

    **Energy axis only.** It is what `positives` reads and what distinguishes
    one lifestyle from another; the carbon axis is an AND over ten separate
    solves and has no single flux solution to count a support over. That is a
    scoping decision, not an omission — see the module docstring.
    """
    from reframed import pFBA

    drain = probe.drains.get("atp")
    substrate = probe.light() if life.key == "photo_co2" else life.substrate
    if drain is None or not substrate or not probe.has(substrate) \
            or not probe.has(life.acceptor):
        return NOT_TESTABLE

    full = tuple(background) + tuple(substrate) + life.acceptor
    without = tuple(background) + life.acceptor
    if probe.maximum(drain, without) > probe.bio.MIN_FLUX:
        return NOT_NEEDED
    ceiling = probe.maximum(drain, full)
    if ceiling <= probe.bio.MIN_FLUX:
        return NO_ROUTE

    # Pinned clear of the solver's tolerance so the route has to carry real
    # flux. `flux_from` reports an unbounded optimum as infinity, and a floor of
    # infinity is an infeasible LP, so that case takes a fixed floor instead.
    constraints = probe.medium(full, drain)
    constraints[drain] = (0.1 * ceiling if math.isfinite(ceiling) else 1.0, None)
    solution = pFBA(probe.model, constraints=constraints, solver=probe.solver)

    gapped = attributable = 0
    for rid, value in (getattr(solution, "values", None) or {}).items():
        if abs(value) <= probe.bio.MIN_FLUX or rid.startswith("R_EX_"):
            continue
        if rid in probe.drains.values() or scores.get(rid, 0.0) > EVIDENCE_MIN:
            continue
        if rid in associable:
            gapped += 1
        else:
            attributable += 1
    return f"{gapped}+{attributable}"


def read(bio, path: Path) -> Probe:
    from reframed import load_cbmodel

    return Probe(bio, load_cbmodel(str(path), flavor="bigg"))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="lifestyle_calibration",
                                description=__doc__.split("\n\n")[0])
    p.add_argument("--workdir", type=Path, required=True,
                   help="where to download and carve into; existing files reused")
    p.add_argument("--organism", action="append", default=[],
                   help=f"restrict to these keys "
                        f"({', '.join(o.key for o in ORGANISMS)}). Repeatable")
    p.add_argument("--curated-only", action="store_true",
                   help="skip carving; score the BiGG models alone")
    p.add_argument("--model", type=Path, action="append", default=[],
                   help="an extra model to score, e.g. a review-set draft. "
                        "Repeatable, and `--organism none --model ...` scores "
                        "nothing else")
    p.add_argument("--faa", type=Path, action="append", default=[],
                   help="an extra proteome to carve and then score. Repeatable")
    p.add_argument("--benchmark", action="store_true",
                   help=f"fetch and score the {len(BENCHMARK)} benchmark genomes "
                        f"({', '.join(b.key for b in BENCHMARK)}) — the set that "
                        "asks whether a column can vary at all")
    p.add_argument("--universe", help="carve from a named CarveMe universe "
                        "(archaea, grampos, gramneg, cyanobacteria) instead of "
                        "the default bacterial one. Overrides the per-organism "
                        "setting, and is `--set carveme--universe=X` in the "
                        "pipeline")
    p.add_argument("--yield-min", type=float, default=YIELD_MIN, metavar="RATIO",
                   help=f"how much more ATP the acceptor must buy before a "
                        f"column calls it respiration (default {YIELD_MIN}, and "
                        "uncalibrated — see Probe.respires)")
    p.add_argument("--evidence", action="store_true",
                   help="also run arm three: restrict the universe to what this "
                        "genome has gene evidence for and count how many "
                        "evidence-free reactions the route still needs")
    args = p.parse_args(argv)

    bio = load_biosynthesis()
    background = salts(bio)
    args.workdir.mkdir(parents=True, exist_ok=True)
    chosen = [o for o in ORGANISMS
              if not args.organism or o.key in args.organism]
    # Empty is legitimate when there are `--model` paths to score: the third
    # question this has to answer is whether a column *varies* across the
    # genomes users actually have, and those models are already carved.
    if not chosen and not args.model and not args.faa and not args.benchmark:
        raise SystemExit(f"lifestyle: no organism matched {args.organism}")

    print(f"salts: {' '.join(background)}\n")
    print("legend: CE substrate gives carbon and energy, E energy only, "
          "C carbon only,\n        - neither, na no exchange to ask through, "
          "? the model cannot be asked\n        that axis at all. gate: ! makes "
          "ATP from nothing, na no ATP drain to test.\n"
          "        evidence rows are k+t — reactions with no gene evidence in "
          "this genome,\n        and reactions no gene could explain. An upper "
          "bound; k=0 is exact.\n")

    header = f"{'model':<26} {'gate':<5} " + " ".join(
        f"{key[:9]:>9}" for key in LIFESTYLE_KEYS)
    print(header)

    rows: dict[tuple[str, str], dict[str, str]] = {}
    gates: dict[tuple[str, str], bool | None] = {}
    voided: dict[tuple[str, str], str] = {}

    def show(label: str, kind: str, gate: str, row: dict[str, str]) -> None:
        print(f"{f'{label}/{kind}':<26} {gate:<5} "
              + " ".join(f"{row[key]:>9}" for key in LIFESTYLE_KEYS))

    def measure(label: str, kind: str, path: Path) -> None:
        probe = read(bio, path)
        leaky = free_energy(probe)
        rows[(label, kind)] = row = score(probe, background, args.yield_min)
        gates[(label, kind)] = leaky
        # What the model could not be asked, named rather than left in the `?`.
        short = []
        if probe.missing_atp:
            short.append("ATP drain needs " + " ".join(probe.missing_atp))
        if probe.missing_precursors:
            short.append("no " + " ".join(probe.missing_precursors))
        if short:
            voided[(label, kind)] = "; ".join(short)
        show(label, kind, "!" if leaky else ("na" if leaky is None else "ok"), row)

    def measure_evidence(label: str, hits: Path, universe_name: str | None) -> None:
        """Arm three: the universe restricted by this genome's own evidence."""
        universe, scores, associable = evidence(hits, universe_name)
        probe = Probe(bio, universe)
        rows[(label, "evidence")] = row = {
            life.key: gap_count(probe, life, background, scores, associable)
            for life in LIFESTYLES}
        show(label, "evidence", "n/a", row)

    for organism in chosen:
        universe_name = args.universe or organism.universe
        sources = [("curated", fetch(BIGG_URL.format(model=organism.bigg),
                                     args.workdir / f"{organism.bigg}.xml",
                                     unzip=True))]
        carved = None
        if not args.curated_only:
            carved = draft(proteome(organism, args.workdir),
                           args.workdir / f"{organism.key}_draft.xml",
                           universe_name)
            sources.append(("draft", carved.model))
        for kind, path in sources:
            measure(organism.key, kind, path)
        if args.evidence and carved is not None:
            measure_evidence(organism.key, carved.hits, universe_name)
    if args.benchmark:
        for entry in BENCHMARK:
            universe_name = args.universe or entry.universe
            carved = draft(benchmark_proteome(entry, args.workdir),
                           args.workdir / f"{entry.key}_bench.xml", universe_name)
            measure(entry.key, "bench", carved.model)
            if args.evidence:
                measure_evidence(entry.key, carved.hits, universe_name)
    for faa in args.faa:
        carved = draft(faa, args.workdir / f"{faa.stem}_draft.xml", args.universe)
        measure(faa.stem, "draft", carved.model)
        if args.evidence:
            measure_evidence(faa.stem, carved.hits, args.universe)
    for path in args.model:
        measure(path.stem, "extra", path)

    print()
    print(f"{'organism':<26} {'curated vs literature':<24} "
          f"{'draft vs curated':<18} notes")
    for organism in chosen:
        curated = rows.get((organism.key, "curated"))
        if curated is None:
            continue
        found = positives(curated)
        hits = found & organism.expected
        literature = (f"{len(hits)}/{len(organism.expected)}"
                      + (f" +{len(found - organism.expected)} extra"
                         if found - organism.expected else ""))
        draft_row = rows.get((organism.key, "draft"))
        if draft_row is None:
            agreement = "—"
        else:
            same = sum(1 for key in LIFESTYLE_KEYS
                       if draft_row[key] == curated[key])
            agreement = f"{same}/{len(LIFESTYLE_KEYS)} cells"
        missed = sorted(organism.expected - found)
        print(f"{organism.name:<26} {literature:<24} {agreement:<18} "
              + (f"curated misses {' '.join(missed)}" if missed else ""))

    if args.benchmark:
        print()
        print(f"{'benchmark genome':<30} {'vs literature':<14} notes")
        for entry in BENCHMARK:
            row = rows.get((entry.key, "bench"))
            if row is None:
                continue
            found = positives(row)
            missed = sorted(entry.expected - found)
            extra = sorted(found - entry.expected)
            print(f"{entry.name:<30} "
                  f"{f'{len(found & entry.expected)}/{len(entry.expected)}':<14} "
                  + (f"misses {' '.join(missed)}" if missed else "")
                  + (f"  claims {' '.join(extra)}" if extra else ""))

    # Gate 4, and the one that can fail after everything else passes: a column
    # that is the same answer for every genome is not a column. That is what
    # took `btn` and `q8` off the biosynthesis panel on 2026-09-10.
    scored = [row for (_, kind), row in rows.items()
              if kind in ("bench", "draft", "extra")]
    if len(scored) > 1:
        print(f"\ndistinct verdicts per column across {len(scored)} drafts — "
              "gate 4 needs two:")
        flat = []
        for key in LIFESTYLE_KEYS:
            seen = sorted({row[key] for row in scored})
            if len(seen) < 2:
                flat.append(key)
            print(f"  {'flat' if len(seen) < 2 else '    '} {key:<11} "
                  f"{len(seen)}  {' '.join(seen)}")
        print(f"  {len(LIFESTYLE_KEYS) - len(flat)} of {len(LIFESTYLE_KEYS)} "
              f"columns vary; flat: {' '.join(flat) if flat else 'none'}")

    leaks = [f"{key}/{kind}" for (key, kind), bad in gates.items() if bad]
    unrun = [f"{key}/{kind}" for (key, kind), bad in gates.items() if bad is None]
    if leaks:
        print(f"\n{len(leaks)} of {len(gates)} models make ATP from nothing — "
              "their energy verdicts are void, not negative:")
        print("  " + " ".join(leaks))
    else:
        print(f"\n{len(gates) - len(unrun)} of {len(gates)} models passed the "
              "free-energy gate.")
    # Separately, because a gate that could not be run is not a gate that passed.
    if unrun:
        print(f"{len(unrun)} could not be gated at all — no ATP drain to "
              "maximise: " + " ".join(unrun))
    if voided:
        print(f"\n{len(voided)} models were short a metabolite a probe needs, so "
              "that axis reads `?` rather than a negative:")
        for (key, kind), why in sorted(voided.items()):
            print(f"  {key}/{kind}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
