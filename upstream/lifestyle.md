# A carbon-and-energy lifestyle table from the metabolic model

Development notes for an idea that is **measured, and half of it is dead**.
Written 2026-09-11, revised 2026-09-15, experiment 1 run 2026-09-15 evening.
The instruments are [`lifestyle_ceiling.py`](lifestyle_ceiling.py) and
[`lifestyle_calibration.py`](lifestyle_calibration.py); this file is why they
exist, what they had to answer, and what a later session must not re-derive
wrongly.

Nothing here is in the pipeline. `grep -rn "lifestyle" src/ docs/ tests/`
returns nothing as of 2026-09-15, and that is now the likely permanent state
for most of it.

## Read this first: it was measured, and the answer is no. Do not rebuild it.

Two experiments ran on 2026-09-15 and between them they close the idea. The
design reasoning below is kept because it is good and a later session will
otherwise re-derive it — but it is building toward a thing that does not work,
and the two reasons are independent, so neither is fixable by patching the
other.

1. **The database cannot support autotrophy, lithotrophy or methanogenesis.**
   Section immediately below.
2. **The acceptor axis, which is what survives, cannot be scored to better than
   a coin flip.** Precision 0.50 on the *curated* models. Section after that.

Full numbers in [../STATUS.md](../STATUS.md); the decision and its reasoning in
[../DECISIONS.md](../DECISIONS.md), the two 2026-09-15 evening entries.

## Experiment 1: the ceiling

**Autotrophy, lithotrophy and methanogenesis cannot be built on CarveMe drafts,
by any scoring, ever.** `CODH_ACS`, `MCR`, `HDR`, `FMFD_b`, `RNF`, `HYDFDN2r`,
`PSI` and `PSII` are in `bigg_universe` and in **none of the five universes
`carve` can load**. A draft is a subnetwork of the universe it came from, so
this is a database limit, and the stopping condition this document set for
itself was exactly that. Arm three does not help: restricting a universe to a
genome's own evidence cannot add a reaction the universe does not contain.

The acceptor axis — O₂, nitrate, fumarate, DMSO, TMAO, sulfate, Fe(III) — is
fully present in the database. That is the surviving half, and it is the half
where the argument against DRAM is weakest, because those are single
well-conserved operons that marker genes already handle. The two halves came out
the wrong way round from what this document predicted.

## Experiment 2: respiration is a yield, and the yield does not transfer

The acceptor axis is what the ceiling left standing. It was scored on the seven
curated models — no carving, so this is the answer key checking itself, and
drafts can only be worse.

**Both binary differentials fail, in opposite directions, on the same organism.**
Removing the donor asks whether the donor is necessary, not whether the acceptor
is used: `iML1515` scores `lac_so4` because lactate is necessary while sulfate
does nothing, so **E. coli is reported a sulfate reducer**. Removing the
acceptor inverts it: E. coli ferments, so no acceptor is ever necessary and it
stops being an aerobe, 2 of 7. E. coli ferments *and* is the textbook organism
for the whole axis.

**The yield ratio fixes it within one organism.** On `iML1515` — ATP with the
acceptor over ATP without, uptake capped at 10.0 — O₂ 2.80, nitrate 3.00, DMSO
and TMAO 1.50, fumarate 1.40, and **sulfate exactly 1.00**. Seven of seven.

**And then it does not transfer.** Over all 31 testable cells in the seven
models, **precision never exceeds 0.50** at any threshold. The cause is
structural: **13 of the 18 expected-false cells return `inf`**, because a donor
that yields nothing alone makes any working acceptor an infinite ratio, and no
finite threshold excludes an infinity. `ac_fe3` fires `inf` for *E. coli*,
*B. subtilis*, *Synechocystis* and *P. putida* — `FE3Ri` is in all five
universes and any model will put electrons on ferric iron offered to it.

So the open question below — "is the ATP yield comparable across models, or only
within? Probably only within" — is answered harder than it was asked. Only
within, and only within organisms that can ferment, which is where the
denominator is finite.

**Everything below this line was written before both results.** It is kept for
the design reasoning, which is still worth having; it is not a plan any more.

## The question

Score a genome's metabolic model on the ways of harvesting carbon and energy —
autotrophy, the terminal acceptors, fermentation, lithotrophy — the way
`biosynthesis.py` scores biosynthetic capability, and put the result in the
report.

**The claim is against DRAM and METABOLIC, not alongside them.** Those tools
map marker genes to KEGG-module completeness. A module reported 80 % complete
says nothing about whether the pathway connects to the rest of metabolism,
whether the cell can regenerate the pathway's electron carriers, or what it
yields. A network can say all three. The intended output is not a checkmark but

> feasible from this medium, **k** reactions in the flux solution have no gene
> evidence in this genome, **N** mol ATP per mol donor

which is module completeness generalised from a linear pathway to a network,
with a yield attached and an explicit count of what had to be assumed. **That
is the contribution**, and the figure that demonstrates it is the measured
disagreement with DRAM: genomes DRAM calls module-complete that are
network-infeasible, and genomes it calls incomplete that close at k = 1.

## Status, 2026-09-15

| | |
| --- | --- |
| `upstream/lifestyle_calibration.py` | arms one to three implemented; **never run** — `reframed` and `carve` are Linux-only on the laptop, and no Linux host was reachable |
| the three experiments below | none has run |
| the ten benchmark genomes | downloaded, `~/postdoc/cm2-lifestyle/genomes/`, 44.0 Mbp and 42,412 proteins over 63 MB. Outside git; `BENCHMARK` in the script holds the accessions, so `--benchmark` re-fetches them from nothing |
| anything in `src/`, `docs/`, `tests/` | untouched, deliberately |
| `DECISIONS.md` | dated entries 2026-09-11 and 2026-09-15 |

### What the 2026-09-15 review changed

- **Two axes that could not be asked said no instead.** A model missing any one
  of the ten precursors could never satisfy the carbon probe, and one missing a
  piece of the ATP drain could never answer the energy probe; both returned `-`.
  That is the conflation `na` exists to prevent, it is silent, and it pointed
  the same way as the prediction on record — so the run would have confirmed
  itself. Now `?`, suffixed, so `E?` reads "energy yes, carbon unaskable". The
  free-energy gate is tri-state for the same reason: a model with no ATP drain
  used to print `ok`.
- **Arm three is implemented** — `evidence()` and `gap_count()`. It was the
  contribution and it was the unwritten arm.
- **`--universe`**, and *M. barkeri* is carved from the archaeal one. Verified:
  `carve.py:280` has `-u/--universe`, and `universe_archaea.xml.gz` is one of
  five CarveMe ships. In the pipeline this is `--set carveme--universe=archaea`
  and it **needs no code change** — `carve_scip.py` forwards unknown arguments
  to `carve` verbatim (`parse_known_args`, line 197).
- **Ten benchmark genomes**, and the donor defect they exposed. See below.

### The donor is baked into the acceptor columns, and that is a defect

Found while assembling the benchmark set, before any solver ran.
*S. oneidensis* MR-1 respires nitrate, fumarate, DMSO, TMAO and Fe(III) — the
entire acceptor axis, the only axis expected to survive gate 4 — and **does not
catabolise glucose**. Every acceptor column pairs its acceptor with glucose, so
MR-1 answers `-` down that axis for a reason that is about the donor.

The shape of the problem: `glc_fum`, `glc_dmso` and `glc_tmao` have exactly one
organism expected to fire them, *E. coli*, and the natural second positive for
all three is the organism the donor excludes. Three of the six columns expected
to ship rest on a single organism.

Fixing it means either standardising on a donor more widely catabolised than
glucose, or crediting an acceptor when any of a small donor panel works. Both
are design decisions and neither is taken.

## What is verified, and what is not

The distinction is load-bearing here, because the argument for the whole idea
rests on one inference that has **not** been checked.

**Verified 2026-09-11**, against live endpoints and the downloaded SBML:

- Eight BiGG models fetch: `iML1515`, `iYO844`, `iJN678`, `iAF692`, `iHN637`,
  `iAF987`, `iJN1463`, `iSynCJ816`. `iMR1_799` and `iSO783` are 404 — there is
  no loadable curated *Shewanella* at those ids.
- Five UniProt proteomes fetch with content. `UP000001656` (*C. ljungdahlii*)
  answers **200 with zero bytes** — superseded — hence `UP000077020` and the
  empty-body guard in `fetch`.
- `iAF692`, `iAF987` and `iHN637` carry no `EX_o2_e` at all; only `iML1515` has
  `EX_dmso_e` and `EX_tmao_e`; `iJN678` spells light `EX_photon_e`.
- 51 demand/sink reactions across the seven curated models, **exactly one with
  a negative lower bound**: `R_SK_pqqA_kt_c` in `iJN1463` at `lb = -1`.

**Verified 2026-09-15**, in the downloaded models:

- `iJN678` (a bacterium) carries `RBPC`, `RBCh`, `PSI`, `PSII`, `PSI_2`.
- `iHN637` (a bacterium) carries `CODH_ACS`, `CODH4`, `FTHFLi`, `MTHFC`,
  `MTHFD`, `MTHFR5`, `RNF`, `HYDFDN2r`, `FDH7` — the whole Wood–Ljungdahl
  pathway with Rnf and a bifurcating hydrogenase.
- `iAF692` (an **archaeon**) carries `MCR`, `HDR`, `FMFD_b`, `CODHr`, `CODH2r`.
- `carveme/<sample>.tsv` is the 12-column DIAMOND-against-BiGG hit table, a
  declared output of the `carveme` rule since it overwrote Bakta's feature
  table (`carve_scip.py`, `catalogue.py`).

**Not verified, and the argument depends on it:** that reactions present in a
BiGG *model* survive into CarveMe's `universe_bacteria`. "The universe is built
from BiGG models" is reasonable and is not a measurement. `universe_ceiling.py`
is exactly the instrument that settles it and **that is why it runs first**.

**Verified 2026-09-15**, against the CarveMe and ReFramed sources on GitHub —
which is the upstream code, not the installed build, and a version check on the
tool environment is still owed:

- `carve` has `-u/--universe` (`carveme/cli/carve.py:280`), and CarveMe ships
  five universes: `bacteria`, `archaea`, `cyanobacteria`, `gramneg`, `grampos`.
- `carve` derives its DIAMOND output from the *input* path — `blast_output =
  os.path.splitext(inputfile)[0] + '.tsv'` — which is why the hits land beside
  the model, and is the same mechanism that once overwrote Bakta's table.
- `reaction_scoring(annotation, gprs, ...)` returns `(scores, gene2gene)` and
  `load_diamond_results(filename)` takes a path, so arm three's inputs are as
  the reproducer in [README.md](README.md) assumes.
- `pFBA(model, objective=None, obj_frac=None, minimize=False, constraints=None,
  reactions=None, solver=None)` — the signature `gap_count` calls.

Still unverified: that `reframed` exposes a loopless FBA; that
`reaction_scoring` actually runs from a saved `hits.tsv` standalone — there is
in-repo precedent but it has not been run against a pipeline output; and that
the **installed** bioconda build ships all five universe files, since CarveMe's
`MANIFEST.in` does not list them and only `universe_bacteria` is exercised
today. One command settles the last one:
`ls $(python -c "import carveme,os;print(os.path.dirname(carveme.__file__))")/data/generated/universe_*`.

## The design

### Two probes per cell, because the axes are orthogonal

A single growth score collapses chemolithoheterotrophy — H₂ for ATP, acetate
for carbon — into a failure of autotrophy. So each cell measures both:

    energy  max flux through `atp_c + h2o_c --> adp_c + pi_c + h_c`
    carbon  how many of ten central precursors the model can produce

Verdicts `CE`, `E`, `C`, `-`, and `na` for a model with no exchange to ask
through. `na` is the panel's `absent` and matters for the same reason: a strict
anaerobe with no `EX_o2_e` was never asked the aerobic question.

### Everything is a difference, never a level

Each probe runs twice — on salts + substrate + acceptor, and again with the
substrate removed — and the substrate is credited only when removing it takes
the answer away. **One mechanism doing three jobs**, which is why it is the
design and not a control bolted on:

1. **The guard against free-energy cycles.** The panel is protected by mass
   balance — a demand on `M_trp__L_c` carries flux only if the atoms came from
   the medium. ATP is not an atom, so a thermodynamically infeasible cycle
   answers an energy question yes for free. A cycle pays the same with the
   donor and without it, so the difference cancels it.
2. **Honesty when the acceptor carries carbon.** Fumarate, DMSO and TMAO are
   all carbon compounds. A plain producibility test on `glc + dmso` credits
   glucose for carbon that may have come from DMSO; with glucose removed the
   precursors are still reachable, so no carbon claim is made.
3. **It makes autotrophy askable.** On `h2 + co2`, dropping H₂ leaves the
   carbon source in place and the precursors unreachable — which is exactly the
   claim "this genome fixes CO₂ and needs the donor to do it".

A MEMOTE ATP-from-nothing gate runs first per model and **voids** the energy
half of a model that fails it rather than scoring it negative.

### Three things that a reasonable first version gets wrong

- **Ten central precursors, not the classic twelve.** `accoa` and `succoa` are
  off the list: a drain on them asks for net CoA synthesis — sulfur,
  pantothenate, a second pathway — not for carbon assimilation. Same trap as
  the `atp_c` note in `biomass_precursors`. The ten that remain are C/H/O/P.
- **Sulfate and ferric iron come out of the background medium.** M9 carries
  `so4` as sulfur source and `fe3` as an iron source, which puts the acceptor
  of a sulfate reducer and an iron reducer into every medium and makes both
  columns untestable. Free once `accoa` is off the list — none of the ten
  contains sulfur — and `fe2` remains as the iron source.
- **A curated model has a second medium.** Closing every `R_EX_*_e` to uptake
  is the whole medium for a CarveMe draft, but BiGG models also ship demand and
  sink reactions, and a sink with a negative lower bound is a compound arriving
  from nowhere. One of 51 is enough. Pinned as sources, left alone as sinks.

### The revision of 2026-09-15: interrogate the evidence, not the draft

The draft model is the wrong object. Carving is a **global parsimony MILP**; a
rare pathway with a handful of moderately-scored reactions loses to that
objective routinely, which is what `btn` and `q8` were — routes the universe
contains and carving does not keep, wrong in 12 of 12 drafts across 8 species.
Asking the draft whether it can do Wood–Ljungdahl asks a question the MILP
already answered on unrelated grounds.

The better object is already on disk. `carveme/<sample>.tsv` holds the DIAMOND
hits for every genome the pipeline has ever processed; `reaction_scoring(hits,
gprs)` turns them into a per-reaction evidence score over the whole universe —
which is what the MILP consumes and then discards. So:

> restrict the universe to reactions this genome has evidence for, ask the
> lifestyle question there, and report how many evidence-free reactions the
> flux solution still needed.

No new tool, no new database, no recompute. The gap count *k* is the score, and
it degrades gracefully where a verdict on a draft collapses to `-`.

## Two positions reversed — do not re-derive them the old way

**Gap-fill distance was dismissed on 2026-09-11 and is rehabilitated.** The
circularity objection — gap-fill to a medium, then test growth on that medium,
and the answer is yes by construction — applies to gap-fill-then-test-binary.
It does not apply to cost-as-score: *k* = 0 and *k* = 12 are different answers
and neither is assumed. The instability objection was about CarveMe's
25,348-reaction global MILP with its degenerate optima (947.4997 against 943;
1,579 reactions against 1,135 on near-identical siblings). A per-lifestyle
gap-fill problem is small, and there degeneracy means "several equally short
routes exist" while the count stays the answer.

**"Fall back to marker genes off bakta" was the wrong recommendation.** It
reimplements what DRAM and METABOLIC already do. The network-level answer is
the point of the exercise, and a marker fallback concedes it. Markers remain a
fair *comparator* — that is the DRAM-disagreement figure — not a fallback.

## Blockers

**The archaeal universe — downgraded on 2026-09-15, not resolved.**
`catalogue.py` still passes no `--universe`, so an archaeon carved by the
pipeline with no flag still comes from the bacterial universe and
methanogenesis is still unreachable by construction. What changed is that the
fix costs nothing: **`--set carveme--universe=archaea` already works**, because
`carve_scip.py` forwards unknown arguments to `carve` verbatim and `carve` has
the flag. No code change, and no dependency on GTDB-Tk's taxonomy — which was
the design cost that made this look expensive. Carl's call, 2026-09-15: that is
the solution for now.

What it does not do is make a **mixed** run correct. `--set` is per-run, not
per-sample, so a set containing both bacteria and archaea gets one universe for
all of them, and the wrong one for some. Selecting per sample is what would
need the taxonomy, and it is deferred, not solved. Worth fixing independently
of the lifestyle idea either way.

**The report's editorial bar is higher than the panel's.** A per-compound
`none` is hedged and local; "not a nitrate reducer" is a claim about the
organism and gets quoted. The panel already needed a release to say when it is
describing the model rather than the genome.

## The experiments, in the order that kills the idea fastest

**1. Universe ceiling on the lifestyle reaction sets.** `universe_ceiling.py`,
against `universe_bacteria` *and* `universe_archaea`, for each of the fifteen
lifestyles. No carving, no downloads beyond the universes, minutes. This is the
one experiment that can kill the idea outright rather than merely disappoint,
and it settles the unverified inference above. If the universe cannot reach
CH₄ from H₂ + CO₂, that is a database limit and no scoring rescues it.

**2. Three-way comparison on the seven organisms.** For each: the curated BiGG
model, the CarveMe draft, and the universe restricted by that genome's own
DIAMOND scores. **If the third arm recovers what the second loses, the thesis
is proven in one table.** The drafts already write their hit files as a side
effect of the carving `lifestyle_calibration.py` does, so the third arm costs
nothing extra to obtain. Implemented 2026-09-15: `--evidence` adds the third
arm, and *k* is an **upper bound** — pFBA minimises total flux, not the count of
evidence-free reactions, so *k* = 0 is exact and *k* = 7 is not. The exact count
is a MILP per cell and is the next thing to build if the bound proves loose.

**3. Does a column vary?** Two sets, and they are different questions:

- `--benchmark`, the ten genomes chosen so that a working probe *must*
  discriminate. This asks whether the probe can discriminate at all.
- `--organism none --model <the review-set drafts>`, the eight from the
  2026-09-11 review set — four *E. faecium*, two *S. aureus*, *E. coli*,
  *M. genitalium*. This asks whether it discriminates on what a user brings.

Running only the second conflates "the probe is blind" with "these genomes are
alike", which on a set of four *E. faecium* is the likelier explanation. A
column that is the same answer for every genome is not a column; that is what
took `btn` and `q8` off the biosynthesis panel. **Experiments 1 and 2 can all
pass and this one still fail.**

## Gates for putting a column in the report

Per column, not in aggregate, because it is per column that it ships. Restated
2026-09-15: the 09-15 revision moved the object of study from the draft to the
evidence-restricted universe and **left these gates written on drafts**, so they
described the arm the revision had just deprecated.

1. The model being scored passes the ATP-from-nothing gate, and could be *asked*
   it — no ATP drain is `na`, not a pass. Otherwise the energy half is void and
   there is nothing to ship.
2. The curated model reproduces its own organism on that column — otherwise
   there is no answer key.
3. On that column, **the evidence-restricted arm matches the curated model for
   all seven** organisms: *k* small where the organism does it, `none` or a
   large *k* where it does not. Draft-versus-curated is still measured and is
   still informative, but it is now the *comparator* — it is what arm three has
   to beat — rather than the gate. If arm three only matches the draft, there
   was no point moving off the draft.
4. **At least two distinct verdicts across the benchmark set, and separately
   across the review set.** The first is necessary — a probe that cannot tell
   *M. barkeri* from *E. coli* is broken. The second is what decides whether the
   column is worth a user's screen.

If the gates pass, the expectation is that **only the acceptor axis ships** —
O₂, nitrate, fumarate, DMSO, TMAO, fermentation. CompareM2's users have
heterotrophs, so the autotrophy and lithotrophy columns are constant-negative
for that audience whether or not the probe works, and fail gate 4 on the
audience rather than on the method.

Shipping touches five files in the usual pattern: a second TSV on the existing
`biosynthesis` spec in `catalogue.py`, a `guidance.py` entry with its numbers
checked against papers, a `report.py` section, `docs/`, and tests.

## Open questions

- Does `reaction_scoring` run standalone off a saved `hits.tsv`, and are the
  scores stable enough to threshold? Nothing has measured what a sensible
  evidence cutoff is.
- What is the right *k*? A route needing 2 transport reactions without evidence
  is different from one needing 2 core enzymes. Partly addressed 2026-09-15:
  `gap_count` reports `k+t`, separating reactions BiGG knows a gene for and this
  genome lacks from reactions no gene could ever explain. What weight the two
  deserve is still unmeasured, and reporting them apart is the least that can be
  done while it is.
- *k* is an upper bound and nothing has measured how loose. If pFBA's support
  routinely carries reactions a minimal solution would not, the number is
  unusable as a score even where the ordering is right — and the fix is the
  per-cell MILP, which is a real cost.
- Is the ATP yield comparable across models with different ETC stoichiometry,
  or only within a model? Probably only within — which would make it a
  within-genome ranking rather than a cross-genome number, and that changes how
  the report can present it.
- No curated sulfate reducer and no curated aerobic methylotroph is in the
  calibration set, so `lac_so4` and `meoh_o2` are pure specificity columns
  *there*: every one of the fourteen models should answer `-` or `na`, and
  anything else is a false positive worth chasing. The benchmark set closes the
  other half — `dvul` and `mext` are the positives those two columns had
  nowhere, and every one of the fifteen columns now has at least one organism
  expected to fire it.
- Which donor the acceptor columns should use. Glucose excludes
  *S. oneidensis*, the one organism in the set with the full acceptor
  repertoire, and three of the six columns expected to ship rest on *E. coli*
  alone as a result. See *the donor is baked into the acceptor columns* above.

## Artifacts

`upstream/lifestyle_calibration.py` — the instrument. Self-contained, runs
under the tool environment's own python, imports `biosynthesis` by path and
nothing from `comparem2`, like `carve_scip.py` and `panel_calibration.py`.
**Activate the environment, do not just call its interpreter**: `carve` shells
out to DIAMOND.

```
conda activate <the carveme env>

# experiment 2, all three arms
python3 upstream/lifestyle_calibration.py --workdir /tmp/lifestyle --evidence

# experiment 3a — can the probe discriminate at all
python3 upstream/lifestyle_calibration.py --workdir /tmp/lifestyle \
    --organism none --benchmark --evidence

# experiment 3b — does it discriminate on what a user brings
python3 upstream/lifestyle_calibration.py --workdir /tmp/lifestyle \
    --organism none --model ~ghrunner/cm2-report-run/*/carveme/*.xml
```

`--benchmark` downloads the ten proteomes itself from the accessions in
`BENCHMARK`, so it needs nothing staged. Copies of both the genomes and the
proteomes are already at `~/postdoc/cm2-lifestyle/genomes/` on the laptop —
outside git, 63 MB — and `--faa <that>/*.faa` uses them instead.

Seven organisms, each a curated model and a proteome. *M. barkeri* is the one
imperfect pair — curated str. Fusaro against a draft of str. MS
(`UP000033033`), because UniProt has no Fusaro proteome, so its draft-versus-
curated count is strain-to-strain as well.

| organism | curated | proteome | expected |
| --- | --- | --- | --- |
| *E. coli* K-12 MG1655 | `iML1515` | CarveMe bundle | aerobic, acetate, nitrate, fumarate, DMSO, TMAO, fermentation |
| *B. subtilis* 168 | `iYO844` | CarveMe bundle | aerobic, acetate, nitrate, fermentation |
| *Synechocystis* PCC 6803 | `iJN678` | `UP000001425` | photoautotrophy, glucose |
| *M. barkeri* | `iAF692` | `UP000033033` | H₂/CO₂, CO, methanol and acetate methanogenesis |
| *C. ljungdahlii* DSM 13528 | `iHN637` | `UP000077020` | H₂/CO₂, CO, sugar fermentation |
| *G. metallireducens* GS-15 | `iAF987` | `UP000007073` | Fe(III) on acetate |
| *P. putida* KT2440 | `iJN1463` | `UP000000556` | aerobic, acetate |

And ten benchmark genomes, `--benchmark`, with no curated model and therefore no
answer key but the literature. Every accession was resolved against the NCBI
datasets API on 2026-09-15 rather than written from memory, and two would have
been wrong that way: *G. metallireducens* GS-15 is `GCF_000012925.1`, and RefSeq
has **no *M. barkeri* Fusaro at all** — the reference is str. MS, which is the
strain this file already carves from UniProt.

| genome | accession | expected |
| --- | --- | --- |
| *E. coli* K-12 MG1655 | `GCF_000005845.2` | aerobic, acetate, nitrate, fumarate, DMSO, TMAO, fermentation |
| *B. subtilis* 168 | `GCF_000009045.1` | aerobic, acetate, nitrate, fermentation |
| *P. putida* KT2440 | `GCF_000007565.2` | aerobic, acetate |
| *Synechocystis* PCC 6803 | `GCF_000009725.1` | photoautotrophy, glucose |
| *C. ljungdahlii* DSM 13528 | `GCF_000143685.1` | H₂/CO₂, CO, fermentation |
| *G. metallireducens* GS-15 | `GCF_000012925.1` | Fe(III) on acetate |
| *M. barkeri* MS | `GCF_000970025.1` | H₂/CO₂, CO, methanol and acetate methanogenesis — **archaeal universe** |
| *N.* (*D.*) *vulgaris* Hildenborough | `GCF_000195755.1` | sulfate on lactate |
| *M. extorquens* AM1 | `GCF_000022685.1` | methanol aerobic, acetate |
| *S. oneidensis* MR-1 | `GCF_000146165.2` | acetate aerobic — and **the acceptor axis it cannot be asked**, above |

44.0 Mbp and 42,412 proteins over the ten, every one a complete genome at 1–5
contigs.

Expectations are properties of the **organism**, from the literature, and
deliberately not of either model — a curated model is free to omit a pathway
its organism has, and a miss there is a fact about the model.
