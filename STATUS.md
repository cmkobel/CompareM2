# CompareM2 v3 — status

What is currently true of a real run. This file changes whenever something is
re-run, which is why it is not in [DESIGN.md](DESIGN.md) — decisions should not
need editing because a tool was verified again.

Last updated **2026-09-03**, from runs on thylakoid.

## Tool verification: 14 of 14

**Verified means executed** on the four *E. faecium* test genomes and producing
correct-looking output. It never means "installed" — see
[DECISIONS.md](DECISIONS.md#a-solve-is-not-a-working-environment).

**All fourteen have now been executed under conda deployment**, which since
2026-09-03 is the only way a tool arrives — the dates below are the first
execution of each command line, mostly under the pixi environment that no longer
carries the tools. See *Two environments, every tool deployed* below for the run
that closed that gap.

| Tool | Verified | Evidence |
| ---- | -------- | -------- |
| seqkit | 09-02 | per-contig lengths and GC; identical stats for the duplicate pair |
| checkm2 | 09-02, 09-03 | 100% complete all four, 0.47% contamination for the duplicate pair; `--database_path` wants the `.dmnd` **file**, not its directory. First run via the since-deleted `--isolated-launcher`; now the one tool with its own deployed environment |
| bakta | 09-02 | 4/4 against db v6.0 light; db `software-min` 1.11 against bakta 1.12.1 |
| amrfinder | 09-02, again 09-03 | 7 / 7 / 11 / 8 genes; database fetched by the pipeline itself. Re-run under `--use-conda` on 09-03 — byte-identical output, and the fetch and analysis rules shared one deployed environment; see below |
| mlst | 09-02 | ST32 / ST32 / ST117 / ST78 — duplicates agree |
| mashtree | 09-02 | duplicate pair at distance 0.00000 |
| treecluster | 09-02 | needed `--threshold`; v2's defaults adopted |
| skani | 09-02 | `-c 70` accepted; duplicates at 100.00% ANI; also writes an `.af` matrix |
| panaroo | 09-02 | 3,780 clusters, 2,091 core; core alignment 1,934,948 bp |
| snp-dists | 09-02 | 0 SNPs between the duplicate pair, identical gap counts |
| fasttree | 09-02 | at `threads=1`; duplicates at 0.0 branch length. `threads=1` is now measured, not assumed — FastTreeMP buys nothing here; see below |
| carveme | 09-02, again 09-04 | 4/4 SBML models. `carve_scip.py` took one *E. faecium* genome from ~9 min and a truncated model to **50 s and 1,743 reactions** — but that fix **does not generalise**: measured 09-04, three of seven *S. aureus* genomes still hit the 600 s ceiling under both SCIP builds, at 1,236–1,263 reactions against 1,609–1,636. See below, twice |
| gtdbtk | 09-02 | all four *E. faecium* at 99.0–99.2% ANI against a 95.0 radius, `ani_screen`; duplicates identical. Six defects had to be fixed first — see below |
| biosynthesis | 09-03 | 13-step Snakemake run, `--until biosynthesis`, exit 0. **3 s a genome**, duplicate pair identical, and byte-identical to the same probe run against the PyPI wheel on macOS. The curated `iML1515` control returns **31 of 32 de novo** — see below |

### GTDB-Tk: six defects, none of which a solve would show
Found 2026-09-02 while making the tool runnable. Each was invisible because the
rule had never been executed, and the first two were described in comments as
though they already existed:

1. **The batchfile was never written.** `--batchfile` named a path no rule
   created, so the command would have died on its first line. Now `Tool.files`.
2. **The summaries were never merged.** GTDB-Tk writes `bac120` and `ar53`
   separately, so the declared `gtdbtk.summary.tsv` could not appear and the
   job would have failed on a missing output even with a working command. Now
   `Tool.post`.
3. **`--skip_ani_screen` does not exist in 2.7.2.** It belonged to the versions
   that screened with Mash; 2.7 screens with skani and dropped the flag.
   `unrecognized arguments` was the only way to find this, short of reading the
   installed tool's `--help` — which is how it was found.

4. **The database was the wrong release**, upstream of all of the above. The
   catalogue pointed at r226; `gtdbtk/config/common.py` in 2.7.2 reads
   `COMPATIBLE_REF_DATA_VERSIONS = ['r232']`, so the tool would have rejected
   the data *after* a 141.4 GB download. r232 is 60.8 GB, less than half,
   because it replaced FastANI's reference genomes with skani sketches.

The last two needed the real run, and both were in the fix for (2):

5. **`python -m comparem2.steps` could not import its own package.** pixi
   provides it through a relative `PYTHONPATH=src`, and a Snakemake rule runs
   with the output directory as its working directory, so `src` did not
   resolve. The classification succeeded and the run then died at the last
   line with `ModuleNotFoundError`. The step is now named by absolute script
   path, and `steps.py` imports nothing from its own package so that it can
   run as a plain script — with a test to keep it that way.
6. **The merge counted the same file twice.** `classify_wf` writes
   `gtdbtk.bac120.summary.tsv` into both `--out_dir` and `--out_dir/classify`.
   Both were passed as candidates because the documentation does not say which
   it uses, and unioning the matches would have merged four genomes into an
   eight-row table. The patterns are candidate *locations* for one file, so
   they now deduplicate by basename.

### GTDB-Tk, end to end
Completed 18:57 on 2026-09-02, `CM2_EXIT=0`, report written. Classification
itself took **17 seconds** for four genomes: r232 ships skani sketches, every
genome cleared the ANI pre-screen, and the identify and align steps were
skipped entirely. The 60.8 GB download was 100% of the wall time.

| | |
| --- | --- |
| all four genomes | `s__Enterococcus_B faecium` — `_B` is GTDB's split genus, not a typo |
| ANI / radius | 99.18, 99.18, 99.13, 99.03 against 95.0 |
| aligned fraction | 0.917, 0.917, 0.908, 0.872 |
| method | `ani_screen` for all four |
| **the duplicate pair** | **identical on every column** — the cross-check holds for the thirteenth tool |

The four questions this file said to ask, answered against the real file:

1. `closest_genome_reference_radius` **exists** under exactly that name, so the
   report's ANI colouring is live rather than silently degraded to plain text.
2. `classification_method` is `ani_screen` — 10 characters, against the
   51-character synthetic worst case. No column-width problem.
3. `user_genome` holds the **batchfile's genome id**, so it joins to every
   other section without translation.
4. The merged file is one header and four rows. Only `bac120` exists for an
   all-bacterial set, and `merge-tsv` tolerated the absent `ar53`.

The summary carries **20 columns, not the 18 the documentation lists** —
`closest_placement_radius` and `translation_table` are extra. Reading columns
by name rather than by position is why that cost nothing.

### The gtdbtk report section is written, and now validated on real output
`_section_gtdbtk` exists and renders: the lineage every genome shares stated
once, then a column per rank where they diverge, ANI coloured against *that
reference's* species radius rather than a global 95%, and a note when a genome
carries no species name. Twelve tests, including a real render.

Rendered from the real summary at 18:57: the shared lineage collapsed to
`Bacteria › Bacillota › Bacilli › Lactobacillales › Enterococcaceae ›
Enterococcus_B`, one Species column, ANI green against the 95.0 radius. The
renderer ran; the generic fallback did not.

It was written before the tool had ever produced a file, so its own tests are
still synthetic; what they could not answer is answered above. The header check
that would have degraded it to the generic table stayed in, and the 20-column
reality it did not expect cost nothing because columns are read by name.

Still untested there: an archaeal genome, so the `ar53` half of the merge and
the two-domain column layout have only ever seen fixtures.

### CarveMe was nine minutes for the wrong reason
Traced 2026-09-02. Nothing in CarveMe's Python is slow. Of the 609.3 s a genome
took (116_2, 2,588 proteins, 24 idle cores): DIAMOND 1.8 s, loading the universe
model 1.7 s, scoring reactions 4.0 s, writing the SBML 0.9 s — and **601 s in
the MILP**, which is CarveMe's own hardcoded 600 s limit expiring. So no rewrite
of the pipeline's or CarveMe's Python could have reached 2% of it.

The MILP is not hard. It is the solver build. Same machine, same interpreter,
same carveme 1.6.6, and a byte-identical problem — written out from each build,
md5 `fba2ad10…` both times:

| SCIP | MILP | status | model | annotated rx dropped | total |
| ---- | ---: | ------ | ----- | -------------------: | ----: |
| conda-forge 10.0.3, as shipped | 601 s | timelimit | 1,193 rx / 750 met | 253 | 609.3 s |
| the same, `presolving/milp/maxrounds=0` | 42.7 s | gaplimit | 1,742 rx / 1,175 met | 44 | 50.8 s |
| PyPI wheel 10.0.2 | 9.8 s | optimal | 1,738 rx / 1,176 met | 45 | 17.6 s |

The wheel does not link PaPILO and the conda-forge builds do
(`printExternalCodeVersions()`). Not the SCIP version: conda-forge 10.0.2
carries PaPILO 3.0.0 and was still in the MILP when it was stopped at 300 s,
seven times the 43 s the same version takes with the presolver off. Not symmetry
handling: `misc/usesymmetry=0` ran to the 900 s limit.

**And the shipped run returns a worse model, not just a later one** — 253 of
1,069 annotated reactions dropped against 45. E8202 (3,185 proteins) repeats
it: 908.1 s and 261 dropped, against 16.7 s and 45.

Reduced to the `scip` command line, no CarveMe and no Python — one build
(conda-forge 10.0.3 / SoPlex 8.0.3 / PaPILO 3.0.1), one written-out problem,
`limits/gap 0.001`:

| run | reported | time | primal | dual |
| --- | -------- | ---: | -----: | ---: |
| `read lp; optimize` | time limit | 300.0 s | 913.500 | **935.696** |
| `read lp; read sol_947; optimize` | gap limit | 4.3 s | **947.500 accepted as feasible** | 947.796 |
| `read lp; read sol_943; optimize` | **optimal** | 5.1 s | 943.500 | 943.500 |
| `read lp; set presolving milp maxrounds 0; optimize` | **optimal** | 7.5 s | **947.500** | 947.500 |

Three different objective values out of one build and one file, two of them
labelled optimal, and a dual bound 11.8 below a point the same build certifies
feasible.

**But "optimal" is not well defined on this problem, so the fix is not framed as
recovering the optimum.** `minmax_reduction` couples each flux to its indicator
with bigM=1e3 against eps=1e-3, and both of those points are feasible only to
tolerance: 5.4e-11 and 5.0e-11 on the constraints, but exactly 1e-6 on
integrality against a default `feastol` of 1e-6 — and **rounding the binaries to
exact integers makes the remaining LP infeasible for both**, checked
independently of the solver's own checker. A SCIP maintainer would be within
their rights to call that a scaling problem rather than a presolve bug, which is
how it should be reported. What survives either reading is the ordering: with
this presolver the run is 12–60x slower and keeps a quarter fewer of the
reactions the annotation supports.

Fixed by `src/comparem2/carve_scip.py`. **Verified through the pipeline on
thylakoid at 19:22**, one genome, `--until carveme`: bakta then carveme, 3/3
steps, exit 0, report written. The carveme rule took **50 s** wall (19:21:21 to
19:22:11), and the report counts **1,743 reactions / 1,175 metabolites / 728
genes** for 116_2. The counts in the table above are of distinct SBML reaction
and species ids and run one reaction lower than the report's own parse.

The log line to grep for when a genome takes nine minutes again is
`carve_scip: presolving/milp/maxrounds=0`.

Worth reporting upstream — every conda install of CarveMe has this.

### And it does not generalise: the slow instances are the genomes, not the build
Measured 2026-09-04 on thylakoid, four solves in parallel on 24 cores, through
`carve_scip.py` exactly as the rule invokes it, against the Bakta proteins the
09-02 seven-genome *S. aureus* run had already produced. `cdcdbc7` now prints
SCIP's own account of each solve, which is where these numbers come from.

conda-forge SCIP 10.0 / pyscipopt 6.2.1, PaPILO linked, presolver off — the
shipped configuration:

| genome | status | gap | MILP | model |
| ------ | ------ | --: | ---: | ----- |
| NCTC8325 | **optimal** | 0 | 10.3 s | 1,609 rx / 1,097 met |
| N315 | timelimit | 1.475% | 600.0 s | 1,246 rx / 833 met |
| MRSA252 | timelimit | 1.430% | 600.0 s | 1,263 rx / 844 met |
| GCF_026920295.1 (draft) | timelimit | 1.482% | 600.0 s | 1,236 rx / 823 met |

PyPI wheel, same pyscipopt and SCIP versions, **no PaPILO** among its external
libraries (SoPlex 8.0.2, CppAD, ZLIB, TinyCThread, AMPL/MP, Nauty, sassy,
Ipopt) — the build the section above proposes as the alternative fix:

| genome | status | gap | MILP | model |
| ------ | ------ | --: | ---: | ----- |
| NCTC8325 | gaplimit | 0.033% | 11.5 s | 1,629 rx / 1,107 met |
| N315 | timelimit | 1.042% | 600.0 s | 1,245 rx / 834 met |
| MRSA252 | timelimit | 2.064% | 600.0 s | 1,251 rx / 828 met |
| GCF_026920295.1 (draft) | timelimit | 1.881% | 600.0 s | 1,242 rx / 828 met |

Five things follow, and four of them close a door:

1. **The A/B is clean.** The conda-forge column reproduces 09-02 exactly — the
   same 1,609 / 1,246 / 1,263 / 1,236 reactions and the same metabolite counts,
   two days later on a fresh invocation. This problem is deterministic.
2. **The wheel does not fix it.** The same three genomes time out, at the same
   ~1,240–1,251 reactions. PaPILO was the whole story on 116_2 and is not the
   story here: **the split is by genome, not by solver build.**
3. **Relaxing `limits/gap` would make it worse.** The gap at cut-off is
   1.0–2.1%, so a 2% tolerance would stop *sooner* on an equally sparse
   incumbent. This looked like the cheap fix before the number existed.
4. **More time is already ruled out** — N315 at 3,600 s (09-02) still returned
   `timelimit`, with 7 reactions fewer than the 600 s run.
5. **The objective barely moves while the model changes by a third.**
   MRSA252's incumbent scores 920.4 against NCTC8325's *converged* 907.1 on the
   wheel — higher — with 1,251 reactions against 1,629. CarveMe's objective does
   not reward reaction count, so its near-optimal region holds networks of very
   different size and which one comes back depends on the search path. That is
   the formulation, not SCIP. Even the converged genome is build-dependent:
   1,609 rx at objective 916.3 under conda-forge, 1,629 at 907.1 under the wheel.

So the reaction deficit is **not** a solver-build artefact, is not bought off
with time, and is not a gap tolerance away. What is left is disclosure — the
status has to reach the report and `biosynthesis`, and does not yet — plus a
question about degenerate optima that belongs upstream.

**Not measured, and it decides how bad this is:** whether the sparse models drop
reactions the annotation supports. The `annotated rx dropped` column above (253
against 44) is what convicted the shipped build on 116_2; the equivalent has
never been computed for these three. Until it is, "1,236 rx against 1,609" says
the models differ, not that the small ones are wrong.

### `presolving/milp/maxrounds` does not exist without PaPILO
Found 2026-09-04 while building the wheel environment for the A/B above, and
fixed in `b6ad30d`. `setParam` raises `KeyError('Not a valid parameter name')`
on a SCIP that does not link PaPILO, so `carve_scip.py` died at the first solve
— exit 1, no model, four for four. The wrapper's own docstring says it "must not
be the reason a run fails", and it had made the build it recommends unusable.
The parameter is now set inside a `try`, and its absence is reported as what it
is: nothing to disable.

### `carve` was overwriting Bakta's feature table
Found 2026-09-02 while tracing the above, and unrelated to speed. `carve`
derives its DIAMOND output path from its *input* path, so
`carve …/bakta/<sample>.faa` wrote `…/bakta/<sample>.tsv`. In the run of that
morning, `samples/116_2/bakta/116_2.tsv` holds 12-column DIAMOND hits against
BiGG gene ids, not Bakta annotations.

Nothing caught it because nothing depends on that file: Bakta declares only its
GFF3 and its FAA, and the report reads the GFF3. The wrapper now links the FAA
into CarveMe's own directory and points `carve` at the link, so the hits land in
`carveme/<sample>.tsv`, where they are a declared output. Confirmed in the 19:22
run: `bakta/116_2.tsv` still opens with `# Annotated with Bakta / # Software:
v1.12.1`, and the DIAMOND hits are in `carveme/116_2.tsv`.

### biosynthesis, executed 2026-09-03 — and what has *not* been executed
`biosynthesis.py` was run directly, against 11 real CarveMe models on this
machine: the four *E. faecium* from `results_full13` and the seven *S. aureus*
from the `verify` run, plus `iML1515` as an external control. 1.2–1.8 s a
genome — negligible against the 50 s `carve` takes — and the `116_2` /
`116_2 duplicate` pair produces byte-identical tables.

| model set | de novo | upstream | no route | absent |
| --- | ---: | ---: | ---: | ---: |
| `iML1515`, curated *E. coli* | **31** | 0 | 1 | 0 |
| 4x *E. faecium* draft | 10–11 | 6–9 | 11–13 | 2 |
| 7x *S. aureus* draft | 25–26 | 0 | 4–5 | 1–3 |

The one non-de-novo verdict on `iML1515` is adenosylcobalamin, which *E. coli*
genuinely cannot synthesise de novo. Every other call on the curated model is
the described answer, which is the calibration for the section.

### biosynthesis through Snakemake, 13:51–14:03 on 2026-09-03
`--until biosynthesis` over the four *E. faecium* genomes: bakta, carveme and
biosynthesis, **13 of 13 steps, exit 0**, report written. The rule takes **3 s**
a genome (14:02:09 → 14:02:12 and the same for the other three) against carveme's
30–70 s, on a machine carrying load 27 from unrelated work.

**The conda build and the PyPI wheel agree exactly.** All four panel tables and
all four media tables are byte-identical to the ones produced on macOS against
`reframed 1.6.0` from PyPI. The conda `reframed` is 1.6.0 too and also reports
only the `scip` backend. The three things that were unverified — that the rule's
`python` is CarveMe's, that conda's ReFramed reads a BiGG-flavoured SBML the same
way, and that `load_cbmodel(..., flavor="bigg")` behaves the same across builds —
are all now answered.

Carving ran fresh rather than reusing the 09-02 models, so this also re-checks
CarveMe: the new `116_2` model has the **same 1,743 reactions, 1,175 metabolites
and 728 genes** as the old one and the same reaction id set, and `E8202` the same
1,813 / 1,221 / 784. The SBML md5s differ and the content does not.

Run from a **fresh clone at `/evo/postdoc/cm2-biosynth-check`**, reusing
`/evo/postdoc/CompareM2/.pixi/envs/default` on PATH rather than solving a second
8.5 GB environment, because the main checkout still carries the pre-commit rsync
state described under *Known broken*.

The 09-03 numbers above were first obtained on **macOS**, which the pipeline
does not support, and only because `reframed` and `pyscipopt` have arm64 wheels
on PyPI where the analysis tools have no `osx-arm64` conda builds. That was a way
to execute this one program's logic early; the Snakemake run is the verification.

### The whole product, produced whole
All thirteen tools in one run over the four *E. faecium* genomes, 19:22 on
2026-09-02: `CM2_EXIT=0`, **95,879 bytes, thirteen sections plus methods**. It
had never happened before — the twelve verified tools were verified across
separate runs and the largest previous report had four sections.

It took two attempts, and the first one is the more useful result. Twelve tools
produced output, amrfinder failed, and **no report was written at all** — see
the two defects below. AMRFinder now reports 7, 7, 11 and 8 genes across the
four genomes, the duplicate pair agreeing.

**339 seconds**, wall, for the whole catalogue over four genomes with every
database already present (19:10:26 to 19:16:05, from the run's Snakemake log —
this is the earlier attempt of the two, which reached carveme and everything
else). Per rule, summed across jobs:

| rule | jobs | total | longest job |
| ---- | ---: | ----: | ----------: |
| bakta | 4 | 350 s | 100 s |
| carveme | 4 | 164 s | 66 s |
| panaroo | 1 | 78 s | 78 s |
| checkm2 | 1 | 52 s | 52 s |
| fasttree | 1 | 27 s | 27 s |
| gtdbtk | 1 | 19 s | 19 s |
| mlst | 1 | 3 s | 3 s |
| seqkit, mashtree, skani, treecluster, snp-dists | 8 | <1 s each | — |

The critical path was bakta → panaroo → fasttree, ending at 339 s; the last
carveme job finished at 300 s, so **carveme is no longer on it**. Before the
presolver fix it would have been the only thing anyone waited for: its four jobs
started at 89, 152, 168 and 234 s and would have run ~610 s each, putting the
end of the run at about 844 s. That figure is arithmetic on this run's observed
start times and the measured 609.3 s per genome — not a re-run.

### Two defects the first full run found
Neither could surface without running the whole catalogue at once.

**AMRFinder's readiness marker outlived its data.** Its database lives in the
tool's conda prefix, not under `--databases`, but the marker recording the
update sat under `--databases` — so rebuilding the environment left a marker
describing data that no longer existed. The download rule was skipped on the
strength of it and the tool failed with `No valid AMRFinder database is found`
against a four-hour-old stamp. An out-of-tree database's marker now lives in
the **run's** output directory, which makes readiness per-run: at most one run
can be wrong about it, and a fresh output directory is always right.

**A partial run produced no report.** `--keep-going` finished twelve tools and
the CLI then returned Snakemake's failure without rendering anything, so the
product was withheld because a thirteenth tool failed. The TUI had always
rendered what succeeded; the CLI now does too, still returning the failure
exit code, and still refusing to write a report for a run that produced
nothing at all.

### Archaea, and both halves of the merge
Two *Methanoflorens* MAGs and two *E. faecium* genomes through gtdbtk, 22
seconds:

| Genome | Domain | Species | ANI |
| --- | --- | --- | --- |
| bog_38 | Archaea | *Methanoflorens stordalenmirensis* | 100.0 |
| fen_3 | Archaea | *Methanoflorens crillii* | 100.0 |
| 116_2 | Bacteria | *Enterococcus_B faecium* | 99.18 |
| E8202 | Bacteria | *Enterococcus_B faecium* | 99.13 |

GTDB-Tk wrote `ar53` and `bac120` summaries separately and `merge-tsv`
combined them into one 5-line, 20-column table — the reason `Tool.post` exists,
finally exercised. The report showed no shared-lineage line and all seven rank
columns, which is the two-domain layout that had only ever seen fixtures.

### The conda deployment model, executed
Both this section and the next describe runs made when `--use-conda` was a flag
and per-rule deployment was one of two models. **The flag was deleted on
2026-09-03 and there are now two environments rather than fourteen** — see *Two
environments, every tool deployed* below. What these two runs established still
holds; the command lines no longer parse.

The thing the whole bioconda package rests on, run for the first time at 19:12
on 2026-09-02 with a **fresh database root**, because the shared one carries an
amrfinder marker a new install would not have.

`--use-conda --until seqkit mashtree treecluster skani checkm2`, exit 0, 10 of
10 steps:

| | |
| --- | --- |
| environments built | **6** — five tools and checkm2's curl+tar download env, as predicted |
| on disk | 5.9 GB in `/evo/postdoc/cm2-conda-envs` |
| addressing | directories named by content hash, e.g. `7af73674e619f38c500a6b917fa7f67a_` |
| **CheckM2** | ran in its own environment **with no `--isolated-launcher`** — the isolation the bioconda package needs, working |
| results | 116_2 and its duplicate both 100.0% complete, 0.47% contaminated; skani 100.00% |
| report | 40,349 bytes, five sections |

That run stopped short of `amrfinder` — the one case the content-addressed
environment sharing actually has to carry. It has now been executed too; see
below.

### amrfinder under `--use-conda`: the environment sharing is real
Run 15:02:33–15:08:40 on 2026-09-03 from `/evo/postdoc/CompareM2` at `e792d5a`,
exit 0, 10 of 10 steps: `download_amrfinder`, 4 × `bakta`, 4 × `amrfinder`.
Four genomes from `tests/E._faecium/`, `--cores 8`, a **fresh `--output`**
(`results_amrfinder_conda`) so the `out_of_tree=True` marker could not
short-circuit the fetch, into the **existing** `--conda-prefix
/evo/postdoc/cm2-conda-envs`. No `--isolated-launcher`.

The claim in `catalogue.py` — that byte-identical env files under one
`--conda-prefix` deploy to one directory, which is what lets the fetch rule and
the analysis rule share a database that can only live inside the environment —
holds, and the log says so rather than the results implying it:

| | |
| --- | --- |
| environments in the prefix | 6 → **8**. Exactly two new: bakta's `e75d6f17a18c7e1351ea2993fbe52345_` and amrfinder's `68e8563502afe8a0983c6c2bb5b459c1_` |
| rules → directories | **4 rules, 2 directories.** `download_amrfinder` and all four `amrfinder` jobs logged `Activating conda environment: ../../cm2-conda-envs/68e8563502afe8a0983c6c2bb5b459c1_` |
| why | the rendered `download_amrfinder.yaml` and `amrfinder.yaml` are byte-identical, md5 `cb5de824e5b0359eeb580a51570bd742` (bakta's pair likewise, `19bdc5c050823a8ef26cbcbd681efcb1`) — the shared `_AMRFINDER_SPEC` doing its job |
| which direction | the **fetch** rule built the environment and the analysis rules reused it, so the sharing does not depend on the analysis rule going first |
| database | 241 MB at `<env>/share/amrfinderplus/data/2026-08-07.1`, `latest` symlinked to it. Version **2026-08-07.1**, as amrfinder prints it — the same version the pixi path fetched |
| fetch cost | **27 s** (25 s download + 2 s `amrfinder_index`), matching the ~26 s estimate |
| amrfinder itself | 7 s per genome, 14 s for four |
| on disk | 8.6 GB in the prefix; amrfinder's environment 1.3 GB of which 241 MB is the database (1010 MB without it), bakta's 1.5 GB |
| marker | `results_amrfinder_conda/amrfinder/.updated`, 0 bytes, 15:08 |
| report | 25,248 bytes, two sections; the AMR section is populated — an SVG heatmap, *10 resistance classes across 4 genomes*, totals 7 / 7 / 11 / 8 |

**The binary that ran came from the deployed environment, not from pixi.** This
needed checking because the pixi environment has amrfinder 4.2.7 *and* database
2026-08-07.1, so identical output would not discriminate. AMRFinder prints its
own paths, and `samples/116_2/logs/amrfinder.log` says
`Software directory: /evo/postdoc/cm2-conda-envs/68e8563502afe8a0983c6c2bb5b459c1_/bin/`
and `Database directory: .../68e8563502afe8a0983c6c2bb5b459c1_/share/amrfinderplus/data/2026-08-07.1`.
Not `/evo/postdoc/CompareM2/.pixi/envs/default/...`, which is what the same
binary reports when run outside an activated environment. The fetch log is the
matching half: `amrfinder_update -d <that same env>/share/amrfinderplus/data`.

Results are **byte-identical to the pixi run** in `results_full13` for all four
genomes — 7, 7, 11, 8 hits. The standing cross-check passes: `116_2` and
`116_2 duplicate` are identical in every column but the protein id, which
carries the sample name.

One operational note that is not about conda-sharing at all: `--use-conda`
needs `conda` on `PATH`, and on thylakoid it is a pixi **global** tool at
`~/.pixi/bin/conda` — nowhere in `.bashrc` or the pixi environment. Without
`export PATH=$HOME/.pixi/bin:$PATH` the run dies at DAG construction with
`Error running conda info. Is conda installed and accessible?`. The 09-02
script had that line; it is load-bearing on this machine and irrelevant to a
real bioconda install, where conda is by definition present.

Script: `/evo/postdoc/cm2-amrfinder-conda.sh`, log
`/evo/postdoc/cm2-amrfinder-conda.log`.

### Two environments, every tool deployed
All fourteen tools with Snakemake installing every one of them, into a **fresh
`--conda-prefix`** so the number is what a new install pays. Run 20:32:27–20:42:01
on 2026-09-03 from `/evo/postdoc/CompareM2` at `b351d3f`, `--cores 16`, exit 0,
**31 of 31 steps**. No `--use-conda` and no `--isolated-launcher`: both had been
deleted, and needing either would have meant the change was wrong.

| | |
| --- | --- |
| environments | **2**, as the CLI predicted before starting (`tool environments: 2`) |
| rules → environments | 30 rule activations: **29 → `main`**, **1 → `checkm2`**. The one is checkm2's own rule |
| on disk | **7.7 GB** — `main` 6.0 GB, `checkm2` 1.8 GB. For contrast the old per-tool prefix is 8.6 GB for **8 of the 14** |
| both solves | **76 s** total (20:32:27 → 20:33:43), main ~49 s and checkm2 ~26 s. **Warm cache** — `~/conda_pkgs_dir` holds 44 GB from earlier runs, so this is not a cold-download figure and should not be quoted as one |
| analysis | 8 min 18 s for four genomes, 14 tools |
| pixi environment | **8.8 GB → 843 MB** once the thirteen tools left `pixi.toml`. It now holds the pipeline, pytest and conda |
| report | 121,912 bytes, **all fourteen tool sections** plus methods |

**The seven tools that had never run this way all did.** Before this, conda
deployment had been executed for seqkit, mashtree, treecluster, skani, checkm2,
bakta and amrfinder. It had not for gtdbtk, mlst, panaroo, snp-dists, fasttree,
carveme or biosynthesis — and two of those were the interesting ones, because
`gtdbtk`'s post step runs our own code through an absolute `sys.executable`
while `carve_scip.py` runs under a bare `python` from the tool's environment.
Both were reasoned about in `catalogue.py` for the conda case without ever
having been executed in it. Both work.

Results agree with the pixi run in `results_full13`, and where they differ the
difference is not in a number:

| | |
| --- | --- |
| byte-identical | seqkit (4/4), amrfinder (4/4), checkm2, gtdbtk, mashtree, treecluster |
| identical but for the embedded absolute path | mlst (ST32/ST32/ST117/ST78), skani (100.00 / 100.00 / 99.14 / 98.91) |
| identical values, different row order | snp-dists — 0 for the duplicate pair, 6,190 to E8202, 8,193 to SRR24, 8,907 E8202↔SRR24. The order comes from panaroo's alignment, which is not deterministic |
| same topology, sixth-decimal branch lengths | fasttree — `((116_2,116_2_duplicate),E8202,SRR24)`, duplicate pair at 0.0, `0.001594939` against `0.001594545` |
| same counts, different cluster names | panaroo — **3,780 clusters, 2,091 core**, matching exactly |
| byte-identical to the 09-03 biosynthesis run | biosynthesis 116_2 — 11 de novo, 6 upstream, 13 none, 2 absent |

The standing cross-check passes on every tool: 116_2 and its duplicate give
2,588 CDS each, 7 AMR genes each, ST32 each, 0 SNPs, 0.0 branch length, 100.0%
complete at 0.47% contamination, identical biosynthesis verdicts, and CarveMe
models 30 bytes apart — 5,602,009 against 5,602,039, which is the sample name
inside the SBML.

Script: `/evo/postdoc/cm2-two-envs.sh`, log `/evo/postdoc/cm2-two-envs.log`.

### `cm2 --setup`, and why installation cannot do it
Added and measured 2026-09-03. Snakemake builds a missing environment during
**DAG construction, before the first job**, so a first run is silent for as long
as the solves take — 76 s on the run above, and on a cluster that is the point
at which someone kills it. `--setup` asks for the same work deliberately.

| | |
| --- | --- |
| cold, fresh prefix | **61.7 s**, both environments, **7.5 GB** |
| against an existing prefix | **1.97 s**, builds nothing |
| assemblies needed | **none** — a 26-byte stub FASTA is written to a temp directory to close the DAG, and deleted afterwards. Nothing reads it because nothing runs |
| databases needed | **none** — `-d` pointed at a directory that does not exist, and it was still not created. Every database path in the DAG is some rule's output, so the graph closes |
| what must match a later run | **`--conda-prefix` only**, however it is set. Its realpath is in the hash; the databases location is not in the env file at all, so `$COMPAREM2_DATABASES` need not be set before `--setup` and may change afterwards. Shown by the acceptance run below, which set up against a nonexistent `-d` and then ran against the real one, reusing the environment |
| tool outputs produced | **zero**, via Snakemake's own `--conda-create-envs-only` |
| scratch left behind | none; `/tmp/comparem2-setup-*` is removed in a `finally` |

**The acceptance test is that a later run reuses it**, and it does: a real
`--until seqkit skani` run against the setup-built prefix activated
`efd5ffa1fbfe3b0c3288c33676aaf20a_` — the same directory `--setup` created —
and finished 4 of 4 steps in **3.5 s** with no environment creation.

**The prefix is in the environment's identity; the output directory is not.**
One byte-identical `main.yaml` deployed to three prefixes on 2026-09-03:

| prefix | directory |
| --- | --- |
| `/evo/postdoc/cm2-envs-two` | `f35bbb1ff167437785dcb4a2729c2beb_` |
| `/evo/tmp/cm2_setup/envs` | `00dbdb482cad56de8230df17a6e17d7c_` |
| `/evo/tmp/cm2_setup_cold/envs` | `efd5ffa1fbfe3b0c3288c33676aaf20a_` |

The other half was measured too: two runs with **different `--output`** and one
prefix, and the second built nothing. That is what lets `--setup` work in a
temp directory it then deletes.

Two footguns that follow. A **relative `--conda-prefix` resolves against
`$INIT_CWD`**, so the same relative path typed from two directories is two
prefixes and two 7.5 GB builds — the shape of the old `./databases` default
bug; the home-relative default is safe, only an explicit relative path bites.
And **`--setup` is one-time per catalogue, not per machine**: bumping any spec
changes `main.yaml`, which changes the hash, which rebuilds all thirteen.

**Which is also why the bioconda recipe cannot do this at install time.** The
prefix is a *runtime* choice — `--conda-prefix`, or `$COMPAREM2_CONDA_PREFIX`,
which on a cluster is set by everyone because home has a quota — so a package
that deployed into `~/.comparem2/envs` during `conda install` would have built
a directory that the hash makes unusable for exactly those users. Two further
reasons, neither measured here and both worth checking before anyone tries:
the only hook is `post-link.sh`, which runs *inside* the conda transaction
installing CompareM2 and would be invoking conda recursively against the same
package cache; and bioconda discourages post-link scripts, expecting them to be
fast, offline-safe and prefix-local. A 7.5 GB, 62 s post-link is none of those.

Where the property *is* available: `cm2 --setup` in a `Dockerfile` layer or a
cluster module's build script. Snakemake also has `--containerize`, which
prints a Dockerfile with the environments baked in — a different thing from the
hand-maintained image rejected on 2026-09-02.

**One cost of having only two environments shows up here.** `Tool.conda` is the
environment's package list rather than the tool's own, so `--until seqkit`
renders the same `main.yaml` as everything: a subset does not get a smaller
environment, it gets the whole 6.0 GB one. `--setup --until <subset>` is
therefore not a cheaper setup. A test records this so it is a known trade rather
than a surprise.

### The standing cross-check
The test set contains `116_2.fna` and `116_2 duplicate.fna` — the same genome
twice, and a filename with a space. **Any tool that treats them differently is
wrong.** Every verified tool agrees: 0.00000 mash distance, 100.00% ANI, 0 SNPs,
identical contig statistics, identical CheckM2 values, model sizes within 30
bytes (the sample name is embedded in the SBML).

### Divergence sanity
A report whose numbers cannot be checked is not a result. 6,190 SNPs between
116_2 and E8202 over a 1,934,948 bp core alignment is 0.320% — the right order
for two strains of one species.

### FastTree's `threads=1` is measured
Checked 2026-09-03, `fasttree 2.2.0` build `h7b50bb2_1` — the build in
`pixi.lock` — with the catalogue's flags, `-nt -gtr`. `bioconda::fasttree` ships
`FastTreeMP` beside `FastTree`, and `Tool.env` can hand it `OMP_NUM_THREADS`, so
the switch was available and was rejected on the numbers.

| input | plain | FastTreeMP | verdict |
| --- | ---: | --- | --- |
| 4 taxa x 1,935,176 bp (test set) | 22.67 s | 23.19 s (x4), 23.23 s (x8) | 2.3–2.5% **slower** |
| 7 taxa x 2,066,459 bp (real run) | 242.3 s, 209.4 s | 251.7 / 278.5 / 254.6 / 204.3 / 208.2 s at 1 / 2 / 4 / 8 / 16 | inside the plain binary's own 14% spread |
| 25 taxa x 100 kbp (simulated) | 25.13 s | 22.31 s (x8) | 1.13x |
| 100 taxa x 100 kbp (simulated) | 132.54 s | 115.67 s (x8) | 1.15x |
| 400 taxa x 100 kbp (simulated) | 597.47 s | 546.18 s (x8) | 1.09x |

The 400-taxon row is two runs of each, ordered plain/MP/MP/plain against a busy
machine. **The payoff does not grow with genome count** — that was the objection
worth testing, and it did not survive. Full numbers and the discarded first pass
are in [DECISIONS.md](DECISIONS.md), 2026-09-03.

Two things that fell out of it, neither acted on:

- **59% of the 7-taxon run is SH-like support** (141.9 s of 239.96 s; 8.5 s of
  21.46 s at 4 taxa), and `report.py`'s `draw_tree` labels leaves only, so those
  values never reach the page. `-nosupport` leaves topology and branch lengths
  byte-identical. At the 4-taxon shape it changes the newick not at all — 3
  unique sequences of 4 leave no non-trivial split to label.
- **FastTreeMP would have shifted branch lengths** in the sixth decimal at 7
  taxa (`0.003849049` → `0.003847998`), from floating-point summation order.
  Topology and support identical, deterministic across thread counts, and at 4
  taxa byte-identical.

## Also verified

- `pixi install` — re-solved 2026-09-03 without the thirteen tools, 8.8 GB
  down to 843 MB, one environment instead of two
- **clean clone** — cloned fresh, environments built, tests pass,
  `test-fast` 8/8 (before the tools left `pixi.toml`)
- Automatic database download — `download_amrfinder` fetched database
  `2026-08-07.1` in 26 s as a workflow step, and again in 27 s from inside a
  deployed environment
- Snakemake lock now lives in the output directory
- **The package builds and installs.** `pip install --no-deps .` into a bare
  venv (Python 3.14, macOS): both entry points present, `comparem2 --version`
  and `cm2 --version` print `CompareM2 3.0.0.dev0`, and the metadata carries
  `License-Expression: GPL-3.0-or-later`.
- **The PATH preflight fired**, on the same installed package with no tools on
  PATH: `not on PATH: seqkit (seqkit), mashtree (mashtree)` and exit 1, before
  any Snakemake call. **That preflight is gone as of 2026-09-03** — nothing is
  expected on PATH any more — and what replaced it is a single check for
  `conda`, which is the failure that actually fires.
- **The deployment flags reach Snakemake unconditionally** — 9.26.1 echoed
  `--software-deployment-method conda --conda-prefix …` and built the DAG. First
  seen under the `--use-conda` flag; now on every run, with no flag.
- **`--tui`, end to end** — `--until mashtree treecluster --tui` under a real
  terminal: two tools selected from the command line, both `done`, report
  written, clean exit. And headlessly through `run_test()` with `seqkit` and
  `checkm2`, so the isolated launcher goes through the TUI path too. Five
  defects had to be fixed first, see [DECISIONS.md](DECISIONS.md). The
  isolated-launcher half of that no longer applies.
- 200 unit tests, ~2.5 s — 8 for the steps around GTDB-Tk's command, 16 for the
  report rewrite, 6 for CarveMe's solver wrapper, 23 for `biosynthesis`, and
  the conda-deployment set rewritten when the flag was deleted
- CI green on GitHub, 18–22 s per run
- `mkdocs build --strict`
- `docs/generate.py --check` — generated pages current

## Environments

**These are the two Snakemake deploys**, measured on `linux-64`. Until
2026-09-03 the same two shapes were a pixi manifest instead; the tool set and
the DIAMOND split are unchanged, which is why the co-solve was not a new risk.

| Environment | Packages | DIAMOND | Contents |
| ----------- | -------: | ------- | -------- |
| `main` | 568 | 2.2.5 | bakta 1.12.1 and twelve other tools, plus curl and tar. 6.0 GB deployed |
| `checkm2` | 127 | 2.1.11 | checkm2 1.1.0 alone. 1.8 GB deployed |

The two DIAMOND versions are the conflict that forced the split, so this is the
isolation working rather than an accident. Both minimum-version pins held —
bakta 1.12.1 and panaroo 1.8.0, not the 1.8.1 and 1.1.2 builds that install
cleanly and crash.

Earlier figure, kept for comparison: 1,142 packages / 1.58 GB, measured with
sylph still in the set.

### macOS on Apple silicon: three packages missing

Solved 2026-09-02 from the laptop, `pixi lock` probes against conda-forge and
bioconda for `osx-arm64` — one per tool plus a combined solve, pixi 0.78.0.
**Nothing was installed and nothing was run**, so this says only what the solver
can reach; it is the weakest kind of evidence in this file. The project stays
Linux-only.

| Missing | Exists for | Required by | Consequence |
| ------- | ---------- | ----------- | ----------- |
| `intbitset` 3.0.2 | conda-forge: linux-64, osx-64, win-64 | panaroo, **every** version 1.5.0–1.8.0 | panaroo uninstallable — the tool is lost, not downgraded |
| `pplacer 1.1.alpha19.*` | bioconda: linux-64 only. arm64 has alpha20 and alpha22 | gtdbtk 2.7.0 / 2.7.1 / 2.7.2 | `catalogue.py`'s `gtdbtk>=2.7` is unsatisfiable; the newest arm build is 2.5.2 |
| `libxcrypt1` 4.4.38 | conda-forge: linux only — it is a glibc crypt library | mlst 2.34.0 and 2.35.0 | mlst falls back to 2.33.1, which does install |

The other ten tools solve. A combined solve of the remaining twelve default-env
dependencies resolves **422 packages** at versions identical to the `linux-64`
lock — seqkit 2.13.0, bakta 1.12.1, amrfinder 4.2.7, mashtree 1.4.6, treecluster
1.0.5, skani 0.3.2, snp-dists 1.2.0, fasttree 2.2.0, carveme 1.6.6 with scip
10.0.3, snakemake 9.16.3, textual 8.2.8. DIAMOND comes out 2.2.6 against linux's
2.2.5. The `checkm2` environment solves too, checkm2 1.1.0 with diamond 2.1.11 —
the same 2.1.x pin that forced the split. These are native arm builds rather
than repackaged x86: build strings and sizes differ from `osx-64`, e.g. skani
637,034 B (`h8f6e10a_0`) against 662,636 B (`h82ec203_0`).

**This probe is what settled the pinning question, and it has since been acted
on.** The solve above used `pixi.toml`'s specs, where `gtdbtk = "*"`; the
combined `osx-arm64` solve did not fail on that, it resolved **gtdbtk 1.0.2
`py_2`, from 2019** — a build that installs cleanly and is junk. That was the
concrete case for floors on every tool rather than three. On 2026-09-03 the
tools left `pixi.toml` entirely and `catalogue.py` became the only place they
are declared, with a `>=` floor on each and a test enforcing it. The probe is
not repeatable as written, because the manifest it probed no longer lists a
tool.

Two things this does not answer: whether gtdbtk 2.5.2 accepts the r232 database
the catalogue downloads (2.7.2 does, by `COMPATIBLE_REF_DATA_VERSIONS`; 2.5.2 is
unchecked), and what CarveMe's MILP costs on arm — every solver number in this
file is from thylakoid.

All three fixes are upstream, none of them local: an `osx-arm64` build on the
intbitset feedstock, a bioconda PR relaxing gtdbtk's stale `alpha19` pin to the
alpha22 that arm already carries, and a report that mlst's *noarch* recipe
declares a dependency on a Linux-only library. **The intbitset one is necessary
but not sufficient** — see below.

### panaroo on arm: intbitset is not the only blocker

Measured 2026-09-03 on the laptop. Unlike the section above, panaroo was built,
installed and **run** here; everything else in this file is from thylakoid.

*intbitset is a packaging gap, not a port.* The source is architecture-neutral —
the only machine-specific construct is `__builtin_popcountll`
(`intbitset_impl.c:157`), and `setup.py` sets `extra_compile_args = []`. Built
from the sdist with conda-forge's arm64 clang: 3.0.2 on py3.11 and 4.1.2 on
py3.13 both compile to `Mach-O 64-bit bundle arm64` and pass **13,930 tests**
(5.93 s, 7.07 s).

The feedstock simply never opted in: `conda-forge.yml` carries no `provider` or
`build_platform` entry for `osx_arm64`, so no rerender will ever emit arm
configs — the bot's newest PR (#20, 4.1.2, 2026-04-28) rerendered and produced
zero. Adding

```yaml
provider:
  osx_arm64: default
```

plus a build-number bump and `conda smithy rerender` emits
`osx_arm64_python3.{11,12,13,14}` — verified against a local clone. The PR has
to carry the version bump too: **3.0.2 fails to compile on py3.12**, while 4.1.2
built as a real `.conda` for arm at py3.11/3.12/3.13/3.14. The changelog from
3.1.0 to 4.1.2 is Python-version support and regenerated C only, no API change,
and panaroo uses the constructor, `|=`, `&`, `len` and iteration.

The feedstock is stalled rather than hostile: last human merge 2023-09-19, seven
open bot PRs, sole maintainer `gtonkinhill` — who also wrote panaroo. This costs
us on Linux today: `python 3.12` + panaroo is unsatisfiable because intbitset
has no py3.12 build, which is what holds the default environment at python
3.11.16.

**The maintainer is active; only the feedstock is stalled.** Checked 2026-09-04,
because the paragraph above reads as though he were gone and he is not: public
activity to 2026-09-03, two panaroo PRs merged 2026-07-02, a push to
`bioconda/bioconda-recipes` on 2026-08-23. That matters for how PR #21 is
handled — it is an argument for waiting rather than escalating.

*The second blocker is prokka, and it cannot be fixed.* With a locally built arm
intbitset in a channel, `panaroo>=1.5` still does not solve: prokka requires
`tbl2asn-forever >=25.7` (linux-64, linux-aarch64, osx-64), which repackages
NCBI's **prebuilt** tbl2asn binary; there is no macOS-arm build of it, and
`table2asn` has none either. prokka 1.15.6 still shells out to it
(`bin/prokka:1429`), and its release notes (2025-12-14) call it the last
release, recommending bakta. Drop prokka and the set solves: **128 packages** on
`osx-arm64`.

panaroo never invokes the prokka binary — `panaroo/prokka.py` is a GFF parser;
only `run_prokka.py:119` shells out, and `run_prokka --help` still exits 0
without it. So the bioconda ask is one deleted line in `requirements: run:`, and
that recipe lists no `recipe-maintainers`, so it is not gated on upstream.

*Executed, not just solved.* panaroo 1.8.0 pip-installed `--no-deps` into that
128-package arm environment; input GFFs from prodigal 2.6.3 over
`tests/E._faecium` (2,587 / 2,587 / 3,192 CDS). `panaroo --clean-mode strict -a
core -t 4 --remove-invalid-genes` exited 0 in **3 min 0 s**: 3,569 clusters,
2,146 core, the duplicate pair differing in **0 of 3,569**. Core alignment
1,952,790 columns; snp-dists 1.2.0 on arm gives **0 SNPs** for the duplicate
pair and 4,111 to E8202. Caveats: prodigal, not bakta, so gene counts are not
comparable to the thylakoid run, and `--remove-invalid-genes` was needed for
prodigal's partial genes (E8202 3,192 → 3,128).

So arm panaroo needs two upstream changes, not one, and neither is a per-arch
recipe: panaroo is `noarch: python` and prokka `noarch: generic`, so
`additional-platforms` is meaningless for both. What is missing is the
dependency closure.

### Both changes are now written, built and proven sufficient

Measured 2026-09-03 on the laptop, later the same day. Notes:
[upstream/intbitset-feedstock-pr.md](upstream/intbitset-feedstock-pr.md) and
[upstream/bioconda-panaroo-pr.md](upstream/bioconda-panaroo-pr.md).

**The intbitset one was sent 2026-09-04 —
[conda-forge/intbitset-feedstock#21](https://github.com/conda-forge/intbitset-feedstock/pull/21).**
The four `osx_arm64` configs built green in conda-forge's own CI on the first
push, which is the claim below confirmed off this laptop. The linter then
failed the PR for a missing `{{ stdlib('c') }}`; adding it needed a second
rerender and a local py3.13 rebuild, and the linter is green after that. The
bioconda one has **not** been sent, and the pair is useless one at a time.

*The cascade is three tools, not one.* Worth stating because the sections above
do not: `snp-dists` and `fasttree` both take
`panaroo/core_gene_alignment.aln` as input (`catalogue.py:512,526`), so
intbitset costs the whole pangenome branch. macOS without these two fixes is
**10 of 14 tools**, not 13.

*The bioconda target moved.* The section above proposed deleting
`tbl2asn-forever` from **prokka's** run deps. That is the wrong recipe: prokka
1.15.6 really does shell out to the binary (`bin/prokka:1429`), so the result
would install and then fail. The change is to **panaroo's** recipe instead —
delete `prokka`. Verified in v1.8.0's source: panaroo's own `setup.py`
`install_requires` never listed prokka, `panaroo/prokka.py` is a GFF parser
with no `subprocess`, and only the separate `run_prokka` console script shells
out. All eight of bioconda's own `test: commands:`, `run_prokka --help`
included, exit 0 with no prokka on PATH.

*The intbitset recipe builds on arm.* `provider: osx_arm64: default` still
emits `osx_arm64_python3.{11,12,13,14}` under conda-smithy **2026.9.1** with
pinning 2026.09.03.11.37.26 — so the finding above holds against current
tooling, and current conda-forge docs confirm `provider`, not `build_platform`,
is the knob. All four configs then build with conda-build 26.7.1, exit 0, real
`.conda` artifacts, peak RSS 228.1 MB. **py3.13 included** — the config the
bot's stale PR #20 fails on. Why #20 fails is *not* known; its logs have
expired (HTTP 410).

*Sufficient, not just necessary.* With those four artifacts in a local channel
and panaroo's recipe carrying the one-line change, the catalogue's own specs
`panaroo>=1.5`, `snp-dists>=1.2.0` and `fasttree>=2.2.0` solve on `osx-arm64`
**from conda alone** — 278 packages, no pip, no `--no-deps`. panaroo 1.8.0
`py_1`, intbitset 4.1.2 `py313h2f2c7d1_0`, **python 3.13.15**, zero prokka.
That python version is the other half of the prize: intbitset 4.1.2 lifts the
cap holding the Linux environment at 3.11.16.

*And it runs.* Installed and executed, because a solve says nothing:

| | pip env (141 pkgs) | conda-only (278 pkgs) |
| --- | --- | --- |
| panaroo wall clock | 180.96 s | 213.27 s |
| gene clusters | 3,569 | 3,569 |
| core (99–100%) | 2,146 | 2,146 |
| shell (15–95%) | 1,423 | 1,423 |
| duplicate pair, clusters differing | **0 of 3,569** | **0 of 3,569** |
| `snp-dists` duplicate pair | **0** (4,111 to E8202) | **0** (4,111 to E8202) |
| `FastTree -nt -gtr` duplicate pair | branch length 0.0 | branch length 0.0 |

`gene_presence_absence.Rtab` is **byte-identical between the two
environments** (sha256 `0d0bebcc…`). `core_gene_alignment.aln` is not, and that
is the record-order nondeterminism panaroo already shows against itself — the
snp-dists distances off both are identical, which is the invariant that
matters. mafft is v7.526 in both, so the 32 s is not the aligner; it is laptop
variation and is not worth a cause.

FastTree had never been run on arm before this. Input is prodigal, not bakta,
so **gene counts still are not comparable to thylakoid** — a bakta-annotated
arm run needs the 1.3 GB light database. That has since been done; see below.

### bakta runs on arm, and gives thylakoid's numbers

Measured 2026-09-04 on the laptop, macOS 26.2 / Apple silicon, 10 cores. Second
tool executed on arm, after panaroo. Working material
`~/postdoc/cm2-macos/bakta/`.

*Installed from conda alone*, `bakta >=1.10` — the catalogue's own spec — into a
pixi project pinned to `osx-arm64`: **155 packages, 1.2 GB**, bakta 1.12.1
`pyhdfd78af_0`, **python 3.13.15**. bakta itself is `noarch`, and every compiled
dependency has a native arm build — `file` reports `Mach-O 64-bit executable
arm64` for diamond 2.2.6, blastn 2.17.0, aragorn 1.2.41, cmscan (infernal
1.1.5), pilercr 1.06 and amrfinder 4.2.7. No pip, no `--no-deps`, nothing built
locally, unlike panaroo.

*Database.* `bakta_db download --type light` — the catalogue's own step —
fetched db **v6.0** (2025-02-24, Zenodo 14916843): 1.34 GB in 4:28 at 5.01 MB/s,
**4.0 GB on disk**, matching the figure in the table below. `version.json`, the
catalogue's `ready` marker, is present.

*Executed*, with `catalogue.py`'s command line verbatim (`--db`, `--threads 8`,
`--output`, `--prefix`, `--force`), on `116_2.fna` and `116_2 duplicate.fna`:
exit 0 both, **66.42 s** and **60.18 s** wall. Not comparable to the 350 s / 100 s
in the per-rule table above — that was four jobs sharing 24 cores, this is two
run one after the other.

| | |
| --- | --- |
| CDS | **2,588** in the GFF3, both genomes — thylakoid's number exactly |
| why the summary says 2,587 | `.txt` counts `CDSs: 2587` and `sORFs: 1` separately; the GFF3 emits both as `CDS`. Same features, two conventions |
| duplicate pair | `.faa` **byte-identical**, sha256 `98c0541b…`; `.gff3` identical but for the lines carrying the sample name. The space in the filename is handled |
| other features | 68 tRNA, 1 tmRNA, 18 rRNA, 13 ncRNA, 34 ncRNA regions, 0 CRISPR, 2 oriC, 305 hypotheticals — identical across the pair |
| AMRFinder expert hits | **7** each, which is what thylakoid reports for these two |

Two caveats. This is bakta alone, not bakta inside a Snakemake run, and no
downstream tool has been fed these annotations on arm — the arm panaroo numbers
above are still prodigal's. And bakta 1.12.1 emits a `SyntaxWarning: invalid
escape sequence '\-'` from `bakta/features/annotation.py:28` under python 3.13;
it is cosmetic, and thylakoid does not see it because intbitset holds that
environment at 3.11.16.

Working material, outside git because it holds a `conda-bld` tree:
`~/postdoc/cm2-macos/` — the prepared feedstock branch, the four `.conda`
artifacts, the applied panaroo patch, and both runs.

## Databases

| Database | Download | On disk | Note |
| -------- | -------: | ------: | ---- |
| GTDB r232 | 60.8 GB | **94 GB** | download measured; on-disk measured 2026-09-03 when the root was moved. 97% of the measured download total. **r226 was wrong**: gtdbtk 2.7.2 accepts only r232 |
| CheckM2 | 1.7 GB | 2.9 GB | measured |
| Bakta light | 1.3 GB | 4.0 GB | download figure is Bakta's documented one; on-disk measured |
| AMRFinder | unmeasured | — | version `2026-08-07.1`, 26 s to fetch. Lands in `$CONDA_PREFIX`, **not** under `--databases` |

Software is 7.7 GB deployed against 62.5 GB of measured downloads, and GTDB-Tk
is 60.8 GB of it. (This line said 1.58 GB until 2026-09-04; that was the
superseded pixi-manifest figure recorded under *Environments* above, measured
with sylph still in the set, and it contradicted the two-environment table on
this same page.) It was 143 GB until the release changed. **On disk the root is
101 GB** (107,812,346,055 bytes, measured when it was moved to `/midifiler`) —
extraction inflates GTDB from 60.8 to 94 GB, which is the figure to plan a
volume around rather than the download size.

## Where things live on thylakoid

24 cores, 125 GB RAM. `/evo` is NVMe, 1.8 T with 785 G free; `/midifiler` is
spinning (`/dev/sda`, `rotational=1`), 13 T with 3.0 T free; `/home` sits on the
other NVMe with 101 G free, which is why neither default location is used.
Measured 2026-09-03. Earlier note: 914 GB free on `/evo` (2026-09-02, before the
60.8 GB download).

| | |
| --- | --- |
| checkout | `/evo/postdoc/CompareM2` — a real clone, current with `origin/master` as of 2026-09-03, and **the only one with a `.pixi`** (843 MB now that the tools are Snakemake's job, down from 8.8 GB). **The path followed the repository rename**; the old lowercase path does not exist |
| scratch clone | `/evo/postdoc/cm2-biosynth-check` — at `e5b181a`, one behind, no `.pixi` of its own; the biosynthesis verification ran here, and it is now deletable |
| databases | `/midifiler/carl/cm2_db_v3` — **101 GB** (94 of it GTDB r232), deliberately **outside** any checkout so deleting a checkout does not cost a re-download. Moved off `/evo` on 2026-09-03 at Carl's call: `/midifiler` is 13 T where `/evo` had 785 G free. 12m23s at 138 MB/s, verified byte-identical at 107,812,346,055 |
| conda prefix | `/evo/postdoc/cm2-envs-two` — 7.7 GB, **2** environments (`main` 6.0 GB, `checkm2` 1.8 GB). The older `cm2-conda-envs` (8.6 GB, 8 single-tool environments) is orphaned and deletable |
| pixi | `/home/thylakoid/.pixi/bin/pixi`. `conda` now comes from the pixi *environment* (a declared dependency, 26.7.1); the pixi **global** at `~/.pixi/bin/conda` is what the 09-03 failure was about and is no longer relied on |

```bash
cd /evo/postdoc/CompareM2

pixi run pytest          # 212 tests, no tools or databases needed
pixi run test-fast       # 4 genomes, no databases needed

pixi run cm2 --setup     # deploy the two environments, once
pixi run cm2 my/*.fna \
  --output results_myrun \
  --cores 24
```

**No path flags needed, as of 2026-09-03.** Both variables are now in
`.bashrc` and a run with neither `-d` nor `--conda-prefix` was verified to pick
them up:

```bash
export COMPAREM2_DATABASES=/midifiler/carl/cm2_db_v3     # 13 T, spinning
export COMPAREM2_CONDA_PREFIX=/evo/postdoc/cm2-envs-two  # NVMe
```

**The split is deliberate.** Databases are 101 GB and mostly GTDB read
sequentially, so they belong on the big spinning volume; the environments are
only 7.5 GB, so moving them would not help, they are tens of thousands of small
files activated by every rule, and moving the prefix **invalidates every
environment** because its realpath is in Snakemake's hash.

**Both exports are load-bearing, and one was actively wrong before this.**
`.bashrc` had three stacked v2-era `COMPAREM2_DATABASES` lines, the last
winning: `/midifiler/carl/comparem2_databases`, which holds 480 GB of *v2* data
and **none of v3's ready markers**. So a `cm2` run without `-d` would have seen
four databases as absent and refetched 62.5 GB into it. Every v3 verification
run passed `-d` explicitly, which is why it never bit. The v2 lines are now
commented; a v2 run needs its path passed. `COMPAREM2_PROFILE` is left alone —
v3 never reads it.

**No `--isolated-launcher` and no `--use-conda`** — both were deleted on
2026-09-03.

To skip the 60.8 GB:

```bash
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme biosynthesis
```

## Packaging

The code side of the bioconda package is done; the release is not. What exists:
`pyproject.toml`, the `comparem2`/`cm2` entry points, `--conda-prefix`, and a
draft recipe in `recipe/`. The model is *pipeline only*, and since 2026-09-03 it
is the *only* model — see [DESIGN.md](DESIGN.md#one-deployment-model-and-two-environments).

| | |
| --- | --- |
| environments a full run builds | **2** — 18 rules, `main` (13 tools + curl + tar) and `checkm2`. Measured at 7.7 GB, built in 76 s warm |
| flags the user needs for any of this | **none.** `--use-conda` and `--isolated-launcher` were deleted. `cm2 --setup` is available to do the build up front |
| published recipe today | `comparem2` **2.16.2**, `noarch: generic`, maintainer `cmkobel` |
| version here | `3.0.0.dev0` — the release needs a real one, in three files |

**The deployment model has been executed whole**: all fourteen tools, both
environments, 31 of 31 steps, correct results, report rendered — see *Two
environments, every tool deployed* above. **There is no untested case left in
it.** The recipe's `run` requirements are already the right shape: the pipeline,
its Snakemake plugins, and `conda`, with the tools deliberately absent.

### The package was built and installed, 2026-09-03

Run in a directory that did not exist beforehand, `/evo/postdoc/cm2-install-test`
on thylakoid, from a `git clone` of `origin/master` at `f9ab013` turned into a
tarball with `git archive --prefix=CompareM2-3.0.0/ HEAD`. That is the same
construction GitHub serves for a tag, so nothing uncommitted in a working tree
could leak into what was tested. **No pixi anywhere in this.**

| step | result |
| --- | --- |
| `conda create` of the recipe's `run` deps + `pip`/`setuptools` | 17.6 s. python 3.12.14, snakemake 9.26.1, textual 8.2.8, setuptools 84.0.0, conda 26.7.1 |
| `python -m pip install . --no-deps --no-build-isolation` — the recipe's build script verbatim | wheel built, `comparem2-3.0.0.dev0-py3-none-any.whl`, 118,589 bytes |
| the recipe's four test commands | all exit 0: `import comparem2`, `comparem2 --version`, `comparem2 --help`, `cm2 --help` |
| `comparem2 --dry-run`, all fourteen tools | **31 of 31 steps** in the DAG, exit 0; the three database rules correctly seen as satisfied |
| `comparem2 --setup` into an empty `--conda-prefix` | **60.9 s**, both environments built. Warm package cache — not a cold-download figure |
| `comparem2 ... --until seqkit mashtree treecluster skani` | 8 of 8 steps, exit 0, **7.4 s**, `report.html` written |
| `conda-build recipe/` | **2m48.9 s**, `comparem2-3.0.0-py_0.conda`, 117,314 bytes, test section passed |
| `conda create ... comparem2` from that package | 16.3 s; `bin/cm2` and `bin/comparem2` both present and working |
| the same four-tool run from the *package* install | 8 of 8, exit 0, **4.9 s**, seqkit output md5-identical to the pip-installed run |

The standing cross-check holds through both installs: mashtree distance 0.00000
and skani 100.00% for the duplicate pair, identical seqkit tables
(`1b7a857b4009c7ee69d822b904e33b4f` from both).

**The one blocker is the version.** The built package is labelled `3.0.0` by
`meta.yaml` while the wheel inside it is `comparem2-3.0.0.dev0.dist-info` and
`comparem2 --version` prints `CompareM2 3.0.0.dev0`. Release step 1 in
[recipe/README.md](recipe/README.md) already covers it; this is what it looks
like if it is skipped.

Four things this answered that a solve could not:

- **`noarch: python` ships no `bin/` files.** `info/files` lists only
  `site-packages/…`; the two entry-point scripts are generated at install time
  from `info/link.json`. Confirmed by installing the built package rather than
  by reading the file list.
- **A read-only package prefix is fine.** `prepare()` writes the Snakefile and
  the env files under `<output>/.comparem2/`, never beside the installed
  module, and both real runs did exactly that. *Superseded 2026-09-04 in one
  respect:* `src/comparem2/` no longer holds only Python. `demo/plasmids.zip`
  is there and does need a `[tool.setuptools.package-data]` entry — measured in
  the built wheel at **572 KB** against 109 KB before, and `--demo` extracts it
  into `<output>/demo_assemblies/` rather than beside the module, so the
  read-only conclusion still holds.
- **`conda` as a run dependency works as intended.** Activating the environment
  puts `conda` on PATH and Snakemake finds it; the deployed environments landed
  in the prefix given by `$COMPAREM2_CONDA_PREFIX`.
- **The recorded `depends` are exactly the five plus python**: `conda`,
  `python >=3.11`, `snakemake-executor-plugin-cluster-generic`,
  `snakemake-executor-plugin-slurm`, `snakemake-minimal >=9,<10`, `textual`.
  License lands as `GPL-3.0-or-later` with `LICENSE` in
  `dist-info/licenses/`.

conda-build emitted two warnings, both benign: the missing source hash (the
source was a local file by design) and its generic "Number of parsed outputs
does not match detected raw metadata blocks", which the two `{% set %}` lines
provoke in recipes that have no `outputs:` at all.

The source tarball is **30 MB, 29 MB of it `tests/`** (the zipped genomes).
Bioconda downloads that on every build. Not a blocker, and not worth splitting
the test data out over.

`/evo/postdoc/cm2-install-test` is **9.8 GB**, 7.7 GB of it a second copy of the
two tool environments in its own prefix. Deletable once this is written down.

### v3.0.0 is on bioconda, and what 3.0.1 owes
Released 2026-09-04. `v3.0.0` points at **`b217b50`**, tarball sha256
**`f7644b3a…`**. [PR #68805](https://github.com/bioconda/bioconda-recipes/pull/68805)
merged at **12:14:21Z** (`8278c92`, approved by bgruening) and
`noarch/comparem2-3.0.0-pyhdfd78af_0.conda` was on anaconda.org at
**12:16:47Z** — 2m26s from merge to available. Tag, merged recipe and published
package all agree; the hash was re-verified against a fresh download and the
61 files of `b217b50` compared byte-for-byte against the tarball.

The tag moved three times — `009c88b` → `b217b50` → `42e2d0d` → `b217b50` —
and the round trip needs explaining, because `b217b50` is knowingly **not** the
best commit. Its own CI is red: `1da59af` added a caveat to `guidance.py`
without regenerating the page derived from it, and the four commits from there
to `ec9b19c` all fail `docs/generate.py --check`. `42e2d0d` fixes that and adds
the unit test that would have caught it.

**It went back because the PR had already merged.** A second recipe update
carrying `42e2d0d`'s hash was pushed to the fork branch at 12:16:51Z — 2m30s
*after* the merge and 4 s after the package was uploaded — so it went nowhere
and the PR head stayed at `ddb0a1c`. (An earlier version of this entry blamed a
GitHub sync failure over eleven minutes. That was wrong: nothing failed to
sync, the branch simply no longer had an open PR.) Moving the tag back to
`b217b50` was then not tidiness but a requirement — the shipped recipe declares
`f7644b3a`, and a tag pointing anywhere else would leave the release
unrebuildable. Re-tagging reproduced `f7644b3a` exactly. The orphan commit
`5e4d7c7` was reset off the branch.

**What went wrong, and it was luck that it did not:** the tag was force-moved to
`42e2d0d` at roughly 12:15–12:16Z, inside the 12:14:21–12:16:47Z window in
which bioconda was building against it. The build succeeded, so it fetched bytes
hashing to `f7644b3a` — but whether that is because it downloaded before the
move or because GitHub served a cached archive is not determinable after the
fact. **Do not move a tag once its recipe is merged and building.** The check
that would have caught it is one line: read `.merged_at` before touching the tag.

The cost of tagging the red commit is bounded and cosmetic: the stale page is
mkdocs source, readthedocs builds from `master`, and the conda package neither
builds nor ships it — three platform builds passing on that exact tarball is the
evidence. What ships is right; what is tagged has a red docs check.

**So 3.0.1 owes the contents of `42e2d0d`** — the regenerated page and its test.
Nothing else is outstanding from this round.

**The published package installs.** `conda install -c conda-forge -c bioconda
comparem2` into a clean environment resolves to v3, checked by Carl on
2026-09-04 — so the `=3` the docs recommend is a guard against a resolver
falling back to a v2 build under conflict, not a workaround for a channel that
does not yet serve v3. An earlier version of the docs said the pin was needed
"until the channel settles", which named no mechanism and was wrong.

Loose end: autobump [PR #68821](https://github.com/bioconda/bioconda-recipes/pull/68821),
opened by the bot at 12:25:44Z after the tag moved. Its hash is correct but it is
a version-only bump — it neither drops v2's `build.sh` nor sets
`noarch: python`, so merging it would undo the recipe #68805 landed. It is on a
branch in the bioconda org, so only bioconda can close it; asked on the PR.

### `cm2 --demo` has been run
Measured on thylakoid 2026-09-04 at `7c98aa9`, against an existing
`--conda-prefix`: **11 of 11 steps, about 3 s** (20:15:08 → 20:15:11), report
**39,716 bytes**. The six bundled plasmids extracted to
`out/demo_assemblies/`, and the seventh input — `116_2 duplicate.fna` —
canonicalised to sample `116_2_duplicate`, space and all.

The pair it exists to demonstrate comes out right in all four analyses:

| analysis | the pair | verdict |
| --- | --- | --- |
| skani | **100.00%** ANI, and identical rows against all five others (93.26 / 95.13 / 92.13 / 93.95 / 94.34) | correct |
| mashtree | **0.00000** branch length | correct |
| treecluster | both in cluster 1 | correct |
| seqkit | `contigs.tsv` md5 `e2ec2407…` for both | correct |

The set is not only a smoke test: `Dallas_55` and `ISMMS_VRE_1` come out at
**99.38%** ANI and share a TreeCluster cluster, while `EF_VRE` sits at 92.13%
against `116_2` — so the ANI matrix and the tree have real structure to read,
which is what makes it worth putting in front of a first-time user.

## Known broken or unfinished

- **AMRFinder's database still lives in the conda prefix**, so it is refetched
  whenever the environment is rebuilt. What is fixed is the *lie*: the marker
  now lives with the run, so a rebuilt environment no longer leaves a stale
  marker claiming the data is there. The refetch was **measured at 27 s** on
  2026-09-03, and the 241 MB it writes sits inside the deployed environment
  rather than under `--databases`. This is a cost, not a defect — see
  *amrfinder under `--use-conda`*. Now shared: the database sits inside `main`,
  which every other tool also uses, so rebuilding it for any reason costs the
  refetch.
- **`--tui` has not been run against a failing workflow interactively.** The
  "Nothing ran / no report" path is covered by unit tests and was reached once
  by accident, but not driven by hand since.
- **`/evo/postdoc/cm2v3`** is the old rsync scratch directory, now redundant,
  holding an 8.5 GB pixi environment that can be deleted.
- **No bioconda package**, and this is now the only thing left: a tag, a
  sha256 and the PR. Every technical unknown in the deployment model has been
  executed — all fourteen tools, two environments, 31 of 31 steps, 2026-09-03.
  See *Packaging* above and `recipe/README.md`. **A hand-built container image
  is no longer planned** — decided 2026-09-02, see
  [DECISIONS.md](DECISIONS.md).
- **The old per-tool conda prefix `/evo/postdoc/cm2-conda-envs` is orphaned.**
  8.6 GB, 8 single-tool environments, addressed by env-file content that no
  longer renders — the two-environment change gives every rule a different hash.
  Deletable. `/evo/postdoc/cm2-envs-two` (7.7 GB) is the live one.
- **The `/evo/postdoc/cm2-databases` copy is still there**, 101 GB, now
  redundant: the live root is `/midifiler/carl/cm2_db_v3`, verified
  byte-identical and exercised by a real checkm2 run and by a run with no path
  flags at all. Deleting it frees 101 GB on the smaller volume, which was the
  point of the move — **not done here because it is 101 GB and GTDB alone is a
  1.7 h refetch**, so it is Carl's call.
- **The thylakoid checkout is current**, fixed on 2026-09-03 before the
  amrfinder run, because that run needed master's `catalogue.py` and the main
  checkout is the only one with a `.pixi`. Carl chose the reversible route over
  the discard the earlier version of this bullet proposed:
  `git stash push -m "pre-v3 rsync snapshot 2026-09-03" -- src tests`, the three
  untracked source files moved to `.rsync-snapshot-backup/`, then
  `git merge --ff-only origin/master`. Now at `12bf984` with only run outputs
  and that backup directory untracked, and the snapshot recoverable with
  `git stash pop`.

  Re-checked against `e792d5a` after the fact, because the earlier version of
  this bullet undersold one file. The tracked nine hold **zero non-blank lines
  master does not have** (`git diff stash@{0} e792d5a -- src tests` is 1,397
  insertions against 19 deletions, and all 19 are blank) — fully redundant. Of
  the three untracked, `steps.py` and `pyproject.toml` are byte-identical to
  master. `carve_scip.py` is **not**, and the claim that it "differs in one
  docstring table header" was wrong: its whole docstring is an earlier
  measurement round — dual bound 934.37 against a feasible 943.4997, and
  `usesymmetry=0` at 907.8 s — where master carries the fuller one, 935.7 at
  300 s against 947.5, the four-configuration comparison, and 253 of 1069
  annotated reactions. Master's is the superset in substance and the more
  careful claim; the backup's prose says the optimum "was presolved away", which
  master retracts as not answerable on this problem. The only precision lost is
  `907.8 s`, which master rounds to the 900 s limit. **Dropping the stash and
  `.rsync-snapshot-backup/` is safe.**
- **`/evo/postdoc/cm2-biosynth-check`** is the clone the biosynthesis
  verification ran from, with `results_biosynth/` inside it. It has no `.pixi`
  of its own — it borrows the main checkout's environment through PATH — and the
  main checkout is now current, so it is **cheap to delete**. It is also one
  commit behind (`e5b181a`), which is why the amrfinder run did not use it.
  Its `results_biosynth/` is still the baseline the 2026-09-03 two-environment
  run's biosynthesis output was compared against, so delete the clone and keep
  nothing else in mind.

## Deliberately left alone

- **The 69 MB of test zips in git.** Real input data, in HEAD, referenced by
  `tests/README.md`. If they ever have to go it is a git-lfs or
  external-hosting decision, not a cleanup one.
- **v2 itself.** Not preserved on a branch, deliberately — a decision taken on
  2026-09-02 when v3 was merged to `master`. `origin/v2` is a 2019-era branch,
  2,008 commits behind the pre-merge `master`, and is *not* what the paper
  describes; the paper-era state is the pre-merge `master` tip, in history and
  on the `ai-1` branch. It is also **tagged**: `origin` carries 39 tags up to
  `v2.16.2` (`ad352f2`), which is what bioconda's published recipe builds from.
  An earlier version of this line said the `v2.*` tags stopped at `v2.9.1` —
  that was this checkout's stale tag list, not the remote's; `git fetch --tags`
  brings the rest. Read the docs at
  [readthedocs `/en/stable/`](https://comparem2.readthedocs.io/en/stable/)
  rather than trying to reconstruct a v2 checkout.
- **The 35.4 MB of PDFs in history.** See
  [DECISIONS.md](DECISIONS.md#lives-on-branch-v3-in-cmkobelcomparem2) — purging
  them would mean rewriting `master`.

### A pixi-installed CompareM2 finds its own conda
Checked 2026-09-04 on macOS, against published 3.0.0 — and worth checking
rather than assuming, because `carve_scip.py`'s neighbouring note records a run
dying at DAG construction with `Error running conda info` when conda was
reachable only through a pixi global.

**There was no conda on `PATH` at all on the test machine** (`command -v conda`
empty), which is what makes the result mean anything. Both pixi routes then ran
`comparem2 <two plasmids> --dry-run --until seqkit` to completion, 3 jobs
planned:

| route | how conda is found |
| --- | --- |
| `pixi add comparem2` in a workspace | `.pixi/envs/default/bin/conda`, on `PATH` inside `pixi run` |
| `pixi global install … comparem2` | `~/.pixi/envs/comparem2/bin/conda`; the trampoline at `~/.pixi/bin/comparem2` sets `path_diff` to that env's `bin`, so the process sees it |

The global install exposes only `comparem2` and `cm2` — **not** `conda` — so the
`path_diff` mechanism is the whole reason it works. If pixi ever stopped setting
it, a global install would break in exactly the documented way.

`pixi global install` needs `--channel conda-forge --channel bioconda`
explicitly unless the channels are already configured; without them it fails
with "No candidates were found for comparem2 *".
