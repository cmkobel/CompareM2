# Upstream reports

Drafts of things this project found in software it depends on, kept here so they
do not evaporate with the session that wrote them. Each file says what it is,
where it goes, and what it still needs. **One has been sent** — the status
column is the record of which.

| draft | destination | status |
| ----- | ----------- | ------ |
| [carveme-205-comment.md](carveme-205-comment.md) | comment on [cdanielmachado/carveme#205](https://github.com/cdanielmachado/carveme/issues/205) | unsent |
| [scip-question.md](scip-question.md) | issue or discussion on [scipopt/scip](https://github.com/scipopt/scip) | unsent |
| [panaroo-intbitset-note.md](panaroo-intbitset-note.md) | issue or discussion on [gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo) | unsent, and incomplete on purpose — see its last section |
| [intbitset-feedstock-pr.md](intbitset-feedstock-pr.md) | PR against [conda-forge/intbitset-feedstock](https://github.com/conda-forge/intbitset-feedstock) | **sent 2026-09-04** — [PR #21](https://github.com/conda-forge/intbitset-feedstock/pull/21) |
| [bioconda-panaroo-pr.md](bioconda-panaroo-pr.md) | PR against [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes) `recipes/panaroo` | unsent; patch prepared and verified to apply |

The first two findings are in [../STATUS.md](../STATUS.md) (*CarveMe was nine
minutes for the wrong reason*), [../DECISIONS.md](../DECISIONS.md) and
`../src/comparem2/carve_scip.py`.

## Six scripts, not drafts

`carve_longsolve.py`, `probe_akg.py` and `find_entry.py` are the instruments
behind the strongest evidence the CarveMe and SCIP drafts have, and neither
draft has been rewritten to use them yet. `panel_calibration.py` and
`universe_ceiling.py` came out of the 2026-09-10 review of `biosynthesis.py`
and are the two checks that should be re-run when CarveMe moves.
`lifestyle_calibration.py` is the newest and the only one that has not been run
at all. All six run under the tool environment's own python and import nothing
from `comparem2`, like `carve_scip.py` and `biosynthesis.py`.

**`lifestyle_ceiling.py` is the seventh, added 2026-09-15, and it is the one
that produced a result.** It asks whether a lifestyle is out of reach for every
possible draft, and the answer is that autotrophy, lithotrophy and
methanogenesis are: `CODH_ACS`, `MCR`, `HDR`, `RNF`, `PSI` and `PSII` are in
`bigg_universe` and in **none of the five universes `carve` can load**. That is
a database limit and no scoring rescues it, which was the stopping condition
[lifestyle.md](lifestyle.md) set for itself. It also found that **the carving
universe is 5,532 reactions, not 25,348** — the larger number is
`bigg_universe.xml.gz`, an input to CarveMe's build that `carve` never loads —
and that `universe_ceiling.py` had been asking that wrong file. 15.8 s to run.

**`panel_calibration.py`** carves the proteomes CarveMe bundles at
`carveme/data/benchmark/fasta/` and scores them, so the panel's calibration is
re-measurable rather than a paragraph: **29 of 30 de novo** for the curated
`iML1515` and for drafts of *E. coli*, *B. subtilis*, *P. aeruginosa*,
*S. oneidensis* and *R. solanacearum* alike, adenosylcobalamin the only miss in
any of them. *M. genitalium* G37 is the negative control at 0 of 30. Needs the
environment activated, not just its interpreter — `carve` shells out to DIAMOND.

**`lifestyle_calibration.py`** has a companion, [lifestyle.md](lifestyle.md) —
the development notes for the idea it measures: the claim against DRAM, what is
verified and what is inference, the two positions reversed on 2026-09-15, the
three experiments in the order that kills the idea fastest, and the gates a
column has to pass to reach the report. Read that first; this paragraph is the
summary.

It asks whether a carved model still knows how its
organism makes a living — the measurement that decides whether a carbon-and-
energy lifestyle table can be built on CarveMe drafts at all. Fourteen models:
seven curated BiGG reconstructions spanning six lifestyles, and drafts of the
same seven organisms. Fifteen lifestyles, two probes each, both read as the
difference between a medium and the same medium with the substrate removed —
which is at once the guard against free-energy cycles, the reason a
carbon-bearing acceptor like DMSO cannot be miscredited, and what makes
autotrophy askable. A MEMOTE ATP-from-nothing gate runs first and voids the
energy half of any model that fails it.

**Written 2026-09-11, arms one to three implemented 2026-09-15, still not run** —
`reframed` and `carve` are Linux-only on this machine, so its Results section is
empty on purpose and its docstring says so. The expectation on record before
running: most lifestyle columns come out `-` for every draft, for the same
reason `btn` and `q8` left the panel — the routes are in the BiGG universe and
carving does not keep them. The nine BiGG and UniProt identifiers it fetches
were checked against the live endpoints on 2026-09-11.

Two things the 09-15 review changed, both of which would have distorted that
run. **An axis that could not be asked said no**: a model short one of the ten
precursors, or missing a piece of the ATP drain, returned `-` rather than a
"could not ask", and the free-energy gate reported `ok` when it had not run.
Sparse drafts are where that bites, and it points the same way as the prediction
— the run would have confirmed itself. And **`--benchmark`**: ten complete RefSeq
genomes, one per way of making a living, accessions resolved against the NCBI
API. They exist because `lac_so4` and `meoh_o2` had no positive anywhere in the
seven-organism set, so those columns could only ever be confirmed negative.

The benchmark set found a defect before it was run. *S. oneidensis* MR-1 has the
whole acceptor axis — the only axis expected to reach the report — and does not
catabolise glucose, which every acceptor column pairs its acceptor with. Details
in [lifestyle.md](lifestyle.md).

**`universe_ceiling.py`** asks whether a panel compound is out of reach for
*every* possible draft, since a draft is a subnetwork of the universe it was
carved from. On 2026-09-10 the answer was that nothing is: all 32 compounds on
the panel as it then stood are reachable, which is what showed `btn` and `q8`
to be carving losses rather than database limits and got them dropped. It is
also the script that exposed `Unbounded` being read as zero flux — the universe
ships every reaction at ±inf, and before the fix this returned 29 of 32 `none`.

**Corrected 2026-09-15: it had been asking `bigg_universe.xml.gz`**, which is
4.6x larger than the universe `carve` loads, so it was not measuring a ceiling.
Re-run against `universe_bacteria` the 09-10 conclusion survives — `btn` and
`q8` are reachable there too — but **`adocbl` is a real ceiling**, absent from
the carving universe entirely. That is why every prototroph misses exactly it
and lands on 29 of 30, curated and carved alike, which nothing had explained.

**`probe_akg.py` and `find_entry.py` found a second, sharper CarveMe defect than
the one `carveme-205-comment.md` currently describes: a draft model that cannot
take up the medium's nitrogen source.** Of the eight *S. mitis* group models
from the 2026-09-08 docs showcase run, **four are missing a link in the
three-step ammonium uptake chain** — three lack `EX_nh4_e` and `NH4tex`
outright, one lacks `NH4tpp` — and those four are exactly the four that report
0 of 32 biosynthesis panel compounds as producible from M9, whose only nitrogen
source is ammonium. The four with a complete chain report 17–18. All eight carry
`GLUDy`, `ASPTA` and `ALATA_L` identically, so the enzymes are not the
difference.

The reason this is a reconstruction bug and not a search accident: the worst
affected model, `Spn_R6`, is the one that reached a **certified optimum in 15.4 s
with the most reactions in the set (1,579)**. CarveMe's objective does not
require the network to be able to eat. This is a cleaner thing to report
upstream than degenerate optima, because it is a yes/no property of the output
that any user can check.

`carve_longsolve.py` is the one that rules out solver time. It re-solves one genome with SCIP's time limit raised, patching
`SCIPSolver.solve` at the same point `carve_scip.py` does. Measured with it on
2026-09-08 (the docs showcase run, *S. mitis* group): **eighteen times the
solver budget bought one reaction.** `Spn_D39` at `limits/time=10800` ran the
full 3 h, still `timelimit` at a 2.01% gap, 1,135 reactions against 1,134 at
600 s — and 444 short of the 1,579 its near-identical sibling `Spn_R6` reaches
in 15.4 s. So the sparse model is not a truncated dense one, which is the claim
`carveme-205-comment.md` currently cannot make. Full numbers in
[../STATUS.md](../STATUS.md), *the docs showcase run*.

## The two panaroo-on-macOS PRs are a pair

`intbitset-feedstock-pr.md` and `bioconda-panaroo-pr.md` are what it takes to
make panaroo installable on Apple silicon, and **neither works alone**. Verified
together on 2026-09-03: with both applied, `panaroo>=1.5`, `snp-dists>=1.2.0`
and `fasttree>=2.2.0` solve on `osx-arm64` from conda alone. Sequence them in
either order, but do not send one and stop.

`panaroo-intbitset-note.md` is a different question — whether panaroo should
carry the dependency at all — and its answer was *probably yes, leave it*. It is
not a prerequisite for either PR.

Working material for both, including the prepared feedstock branch, the applied
patch and the local arm builds, is in `~/postdoc/cm2-macos/` — outside git,
because it holds a conda-bld tree.

Deliberately **not** filed: conda-forge's `scip` feedstock and bioconda's
`carveme` recipe. PaPILO is a legitimate part of the feedstock's build, there is
no PaPILO-free conda SCIP to pin to, and the fix belongs in CarveMe.

## The reproducer

Both drafts offer three attachments: the written-out MILP and two solutions.
They are generated files, not kept in git. To rebuild them from a protein FASTA
in an environment that has carveme:

```python
# write the MILP exactly as `carve` would build it, then stop
import pandas as pd
from carveme import config, project_dir
from carveme.reconstruction.carving import minmax_reduction
from carveme.reconstruction.diamond import load_diamond_results
from carveme.reconstruction.scoring import reaction_scoring
from reframed import load_cbmodel
from reframed.core.transformation import apply_bounds
from reframed.solvers.scip_solver import SCIPSolver

universe = load_cbmodel(project_dir + config.get('generated', 'default_universe'),
                        flavor='bigg')
universe.id = "probe"
apply_bounds(universe)
gprs = pd.read_csv(project_dir + config.get('generated', 'bigg_gprs'))
gprs = gprs[gprs.reaction.isin(list(universe.reactions))]
scores, _ = reaction_scoring(load_diamond_results("hits.tsv"), gprs)

def solve(self, *a, **k):
    self.problem.writeProblem("carveme.lp")
    raise SystemExit(0)

SCIPSolver.solve = solve
minmax_reduction(universe, dict(scores[['reaction', 'normalized_score']].values),
                 default_score=-1.0, uptake_score=0.0, soft_score=1.0)
```

`hits.tsv` is what `carve` leaves beside its input: `diamond blastp -d
<carveme>/data/generated/bigg_proteins.dmnd -q proteins.faa -o hits.tsv
--more-sensitive --top 10 --quiet`.

The two solutions come from the LP itself, no CarveMe needed:

```bash
scip -c "read carveme.lp set presolving milp maxrounds 0 set limits gap 0.001 \
         optimize write solution carveme_947.sol quit"     # the 947.4997 point
scip -c "read carveme.lp set limits gap 0.001 optimize \
         write solution carveme_943.sol quit"              # with the wheel's SCIP
```

The copies used for the numbers in these drafts were produced on thylakoid from
`tests/E._faecium/116_2.fna` on 2026-09-02 and staged, gzipped, in
`~/Downloads/carveme-scip-repro/` on the laptop — which is not a durable
location, hence the recipe above.
