# CompareM2 v3 — status

What is currently true of a real run. This file changes whenever something is
re-run, which is why it is not in [DESIGN.md](DESIGN.md) — decisions should not
need editing because a tool was verified again.

Last updated **2026-09-03**, from runs on thylakoid.

## Tool verification: 14 of 14

**Verified means executed** on the four *E. faecium* test genomes and producing
correct-looking output. It never means "installed" — see
[DECISIONS.md](DECISIONS.md#a-solve-is-not-a-working-environment).

| Tool | Verified | Evidence |
| ---- | -------- | -------- |
| seqkit | 09-02 | per-contig lengths and GC; identical stats for the duplicate pair |
| checkm2 | 09-02 | 100% complete all four, via `--isolated-launcher`; `--database_path` wants the `.dmnd` **file**, not its directory |
| bakta | 09-02 | 4/4 against db v6.0 light; db `software-min` 1.11 against bakta 1.12.1 |
| amrfinder | 09-02, again 09-03 | 7 / 7 / 11 / 8 genes; database fetched by the pipeline itself. Re-run under `--use-conda` on 09-03 — byte-identical output, and the fetch and analysis rules shared one deployed environment; see below |
| mlst | 09-02 | ST32 / ST32 / ST117 / ST78 — duplicates agree |
| mashtree | 09-02 | duplicate pair at distance 0.00000 |
| treecluster | 09-02 | needed `--threshold`; v2's defaults adopted |
| skani | 09-02 | `-c 70` accepted; duplicates at 100.00% ANI; also writes an `.af` matrix |
| panaroo | 09-02 | 3,780 clusters, 2,091 core; core alignment 1,934,948 bp |
| snp-dists | 09-02 | 0 SNPs between the duplicate pair, identical gap counts |
| fasttree | 09-02 | at `threads=1`; duplicates at 0.0 branch length. `threads=1` is now measured, not assumed — FastTreeMP buys nothing here; see below |
| carveme | 09-02 | 4/4 SBML models — but the ~9 min per genome was a solver defect and the models were truncated. **50 s and 1,743 reactions** through `carve_scip.py`; see below |
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

- `pixi install` — both environments, from a manifest that had never been
  installed
- **clean clone** — cloned fresh, both environments built, tests pass,
  `test-fast` 8/8
- `--isolated-launcher` — CheckM2 from its own environment
- Automatic database download — `download_amrfinder` fetched database
  `2026-08-07.1` in 26 s as a workflow step
- Snakemake lock now lives in the output directory
- **The package builds and installs.** `pip install --no-deps .` into a bare
  venv (Python 3.14, macOS): both entry points present, `comparem2 --version`
  and `cm2 --version` print `CompareM2 3.0.0.dev0`, and the metadata carries
  `License-Expression: GPL-3.0-or-later`.
- **The PATH preflight fires.** The same installed package, no tools on PATH:
  `not on PATH: seqkit (seqkit), mashtree (mashtree)` and exit 1, before any
  Snakemake call.
- **`--use-conda` reaches Snakemake with the right flags** — Snakemake 9.26.1
  echoed `--software-deployment-method conda --conda-prefix …` and began
  building the DAG. See below for where that run stops.
- **`--tui`, end to end** — `--until mashtree treecluster --tui` under a real
  terminal: two tools selected from the command line, both `done`, report
  written, clean exit. And headlessly through `run_test()` with `seqkit` and
  `checkm2`, so the isolated launcher goes through the TUI path too. Five
  defects had to be fixed first, see [DECISIONS.md](DECISIONS.md).
- 196 unit tests, ~1.3 s — 13 added for the conda-deployment path, 8 for the
  steps around GTDB-Tk's command, 16 for the report rewrite, 6 for CarveMe's
  solver wrapper, 23 for `biosynthesis`
- CI green on GitHub, 18–22 s per run
- `mkdocs build --strict`
- `docs/generate.py --check` — generated pages current

## Environments

Re-solved after sylph was removed. Measured, `linux-64`:

| Environment | Packages | DIAMOND | Contents |
| ----------- | -------: | ------- | -------- |
| `default` | 568 | 2.2.5 | bakta 1.12.1 and eleven other tools |
| `checkm2` | 127 | 2.1.11 | checkm2 1.1.0 alone |

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

**The pins would have to move into `pixi.toml` first.** With `gtdbtk = "*"` as
the manifest spells it today, the combined `osx-arm64` solve does not fail — it
resolves **gtdbtk 1.0.2 `py_2`, from 2019**. That is the failure mode the
manifest's own comment warns about, and on macOS it would install cleanly and be
junk.

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

## Databases

| Database | Download | On disk | Note |
| -------- | -------: | ------: | ---- |
| GTDB r232 | 60.8 GB | — | measured; 97% of the measured total. **r226 was wrong**: gtdbtk 2.7.2 accepts only r232 |
| CheckM2 | 1.7 GB | 2.9 GB | measured |
| Bakta light | 1.3 GB | 4.0 GB | download figure is Bakta's documented one; on-disk measured |
| AMRFinder | unmeasured | — | version `2026-08-07.1`, 26 s to fetch. Lands in `$CONDA_PREFIX`, **not** under `--databases` |

Software is 1.58 GB against 62.5 GB of measured data, and GTDB-Tk is 60.8 GB
of it. It was 143 GB until the release changed.

## Where things live on thylakoid

24 cores, 125 GB RAM, 914 GB free on `/evo` (measured 2026-09-02, before the
60.8 GB download).

| | |
| --- | --- |
| checkout | `/evo/postdoc/CompareM2` — a real clone, at `e792d5a` = `origin/master` as of 2026-09-03, and **the only one with a `.pixi`**. **The path followed the repository rename**; the old lowercase path does not exist |
| scratch clone | `/evo/postdoc/cm2-biosynth-check` — at `e5b181a`, one behind, no `.pixi` of its own; the biosynthesis verification ran here, and it is now deletable |
| databases | `/evo/postdoc/cm2-databases` — 6.9 GB, deliberately **outside** any checkout so deleting a checkout does not cost a re-download |
| conda prefix | `/evo/postdoc/cm2-conda-envs` — 8.6 GB, 8 environments, shared across `--use-conda` runs |
| pixi | `/home/thylakoid/.pixi/bin/pixi`, and **`conda` is a pixi global** at `/home/thylakoid/.pixi/bin/conda` — `--use-conda` needs that directory on `PATH` |

```bash
cd /evo/postdoc/CompareM2
export PATH=$HOME/.pixi/bin:$PATH
export COMPAREM2_DATABASES=/evo/postdoc/cm2-databases

pixi run pytest          # 196 tests, no tools or databases needed
pixi run test-fast       # 4 genomes, no databases needed

pixi run cm2 my/*.fna \
  --output results_myrun \
  --cores 24 \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

**The export belongs in `.bashrc`, and the launcher path must be absolute.**
The default database directory is `~/.comparem2/databases`, which on thylakoid
is under `/home` and not the `/evo` volume the existing copy sits on, so
without `$COMPAREM2_DATABASES` (or `-d`, which has to be retyped every run) the
6.9 GB is fetched a second time, into home. Without an **absolute**
launcher path, CheckM2 fails with `command not found`, because Snakemake's
shell does not inherit an interactive PATH.

To skip the 60.8 GB:

```bash
--until seqkit checkm2 bakta amrfinder mlst mashtree treecluster skani \
        panaroo snp-dists fasttree carveme biosynthesis
```

## Packaging

The code side of the bioconda package is done; the release is not. What exists:
`pyproject.toml`, the `comparem2`/`cm2` entry points, `--use-conda` and
`--conda-prefix` through both execution paths, and a draft recipe in
`recipe/`. The model is *pipeline only* — see
[DESIGN.md](DESIGN.md#two-deployment-models).

| | |
| --- | --- |
| environments a full conda run builds | **14** (18 rules needing one, deduplicated by content — `biosynthesis` shares CarveMe's) |
| published recipe today | `comparem2` **2.16.2**, `noarch: generic`, maintainer `cmkobel` |
| version here | `3.0.0.dev0` — the release needs a real one, in three files |

**A conda deployment has now been executed** — see *The conda deployment
model, executed* above. Six environments for a five-tool subset, CheckM2
isolated without a launcher, correct results, report rendered. **And
`amrfinder` with it**, the one tool whose database lands inside its deployed
environment: the fetch rule and the four analysis rules share one
content-addressed directory, results byte-identical to the pixi run — see
*amrfinder under `--use-conda`* above. **There is no untested case left in the
conda deployment model.**

## Known broken or unfinished

- **AMRFinder's database still lives in the conda prefix**, so it is refetched
  whenever the environment is rebuilt. What is fixed is the *lie*: the marker
  now lives with the run, so a rebuilt environment no longer leaves a stale
  marker claiming the data is there. The refetch was **measured at 27 s** on
  2026-09-03, and the 241 MB it writes sits inside the deployed environment
  rather than under `--databases`. This is a cost, not a defect — see
  *amrfinder under `--use-conda`*.
- **`--tui` has not been run against a failing workflow interactively.** The
  "Nothing ran / no report" path is covered by unit tests and was reached once
  by accident, but not driven by hand since.
- **`/evo/postdoc/cm2v3`** is the old rsync scratch directory, now redundant,
  holding an 8.5 GB pixi environment that can be deleted.
- **No bioconda package**, but the blocker has moved again: what remains is a
  tag, a sha256 and the PR. The `--use-conda` run that was the last technical
  unknown has been done. See *Packaging* above and `recipe/README.md`. **A hand-built container image is no longer planned** —
  decided 2026-09-02, see [DECISIONS.md](DECISIONS.md).
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
