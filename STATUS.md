# CompareM2 v3 — status

What is currently true of a real run. This file changes whenever something is
re-run, which is why it is not in [DESIGN.md](DESIGN.md) — decisions should not
need editing because a tool was verified again.

Last updated **2026-09-02**, from runs on thylakoid.

## Tool verification: 13 of 13

**Verified means executed** on the four *E. faecium* test genomes and producing
correct-looking output. It never means "installed" — see
[DECISIONS.md](DECISIONS.md#a-solve-is-not-a-working-environment).

| Tool | Verified | Evidence |
| ---- | -------- | -------- |
| seqkit | 09-02 | per-contig lengths and GC; identical stats for the duplicate pair |
| checkm2 | 09-02 | 100% complete all four, via `--isolated-launcher`; `--database_path` wants the `.dmnd` **file**, not its directory |
| bakta | 09-02 | 4/4 against db v6.0 light; db `software-min` 1.11 against bakta 1.12.1 |
| amrfinder | 09-02 | 7 and 11 genes; database fetched by the pipeline itself |
| mlst | 09-02 | ST32 / ST32 / ST117 / ST78 — duplicates agree |
| mashtree | 09-02 | duplicate pair at distance 0.00000 |
| treecluster | 09-02 | needed `--threshold`; v2's defaults adopted |
| skani | 09-02 | `-c 70` accepted; duplicates at 100.00% ANI; also writes an `.af` matrix |
| panaroo | 09-02 | 3,780 clusters, 2,091 core; core alignment 1,934,948 bp |
| snp-dists | 09-02 | 0 SNPs between the duplicate pair, identical gap counts |
| fasttree | 09-02 | at `threads=1`; duplicates at 0.0 branch length |
| carveme | 09-02 | 4/4 SBML models — but the ~9 min per genome was a solver defect and the models were truncated. **50 s and 1,743 reactions** through `carve_scip.py`; see below |
| gtdbtk | 09-02 | all four *E. faecium* at 99.0–99.2% ANI against a 95.0 radius, `ani_screen`; duplicates identical. Six defects had to be fixed first — see below |

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
carries PaPILO 3.0.0 and never finished either. Not symmetry handling:
`misc/usesymmetry=0` still hit the limit, at 907.8 s.

**PaPILO is wrong here, not slow.** The shipped run's dual bound descends to
934.37 while a point of objective 943.4997 exists — and that same build's own
feasibility checker accepts the point, whose largest constraint violation is
4.99e-11 against a `feastol` of 1e-6. An upper bound below a feasible objective
means the optimum was presolved away, which is why the fast path is also the
better model rather than a trade. The likely trigger is CarveMe's conditioning:
`minmax_reduction` couples each flux to its indicator with bigM=1e3 against
eps=1e-3. E8202 (3,185 proteins) repeats all of it: 908.1 s and 261 annotated
reactions dropped, against 16.7 s and 45.

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

### The whole product, produced whole
All thirteen tools in one run over the four *E. faecium* genomes, 19:22 on
2026-09-02: `CM2_EXIT=0`, **95,879 bytes, thirteen sections plus methods**. It
had never happened before — the twelve verified tools were verified across
separate runs and the largest previous report had four sections.

It took two attempts, and the first one is the more useful result. Twelve tools
produced output, amrfinder failed, and **no report was written at all** — see
the two defects below. AMRFinder now reports 7, 7, 11 and 8 genes across the
four genomes, the duplicate pair agreeing.

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

Still untested there: `amrfinder` under `--use-conda`, which is the case the
content-addressed environment sharing in DESIGN.md is reasoned about but not
yet demonstrated.

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
- 167 unit tests, ~1.3 s — 13 added for the conda-deployment path, 8 for the
  steps around GTDB-Tk's command, 16 for the report rewrite, 6 for CarveMe's
  solver wrapper
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
| checkout | `/evo/postdoc/CompareM2` — a real clone, on `master`. **The path followed the repository rename**; the old lowercase path does not exist |
| databases | `/evo/postdoc/cm2-databases` — 6.9 GB, deliberately **outside** any checkout so deleting a checkout does not cost a re-download |
| pixi | `/home/thylakoid/.pixi/bin/pixi` |

```bash
cd /evo/postdoc/CompareM2
export PATH=$HOME/.pixi/bin:$PATH
export COMPAREM2_DATABASES=/evo/postdoc/cm2-databases

pixi run pytest          # 167 tests, no tools or databases needed
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
        panaroo snp-dists fasttree carveme
```

## Packaging

The code side of the bioconda package is done; the release is not. What exists:
`pyproject.toml`, the `comparem2`/`cm2` entry points, `--use-conda` and
`--conda-prefix` through both execution paths, and a draft recipe in
`recipe/`. The model is *pipeline only* — see
[DESIGN.md](DESIGN.md#two-deployment-models).

| | |
| --- | --- |
| environments a full conda run builds | **14** (17 rules needing one, deduplicated by content) |
| published recipe today | `comparem2` **2.16.2**, `noarch: generic`, maintainer `cmkobel` |
| version here | `3.0.0.dev0` — the release needs a real one, in three files |

**A conda deployment has now been executed** — see *The conda deployment
model, executed* above. Six environments for a five-tool subset, CheckM2
isolated without a launcher, correct results, report rendered. What remains
untested there is `amrfinder`, whose database lands inside its deployed
environment; the content-addressed sharing that should make that work is
reasoned about in DESIGN.md and not yet demonstrated.

## Known broken or unfinished

- **AMRFinder's database still lives in the conda prefix**, so it is refetched
  whenever the environment is rebuilt. What is fixed is the *lie*: the marker
  now lives with the run, so a rebuilt environment no longer leaves a stale
  marker claiming the data is there. The refetch costs about 26 s.
- **`amrfinder` has not been run under `--use-conda`.** It is the one tool
  whose database lands inside its deployed environment, which is what the
  content-addressed environment sharing in DESIGN.md is for.
- **`--tui` has not been run against a failing workflow interactively.** The
  "Nothing ran / no report" path is covered by unit tests and was reached once
  by accident, but not driven by hand since.
- **`/evo/postdoc/cm2v3`** is the old rsync scratch directory, now redundant,
  holding an 8.5 GB pixi environment that can be deleted.
- **No bioconda package**, but the blocker has moved again: what remains is a
  tag, a sha256 and the PR. The `--use-conda` run that was the last technical
  unknown has been done. See *Packaging* above and `recipe/README.md`. **A hand-built container image is no longer planned** —
  decided 2026-09-02, see [DECISIONS.md](DECISIONS.md).
- **The thylakoid checkout carries uncommitted `src/` changes**, rsynced from
  the laptop for the GTDB-Tk run rather than pulled. Once the work is committed,
  `git checkout . && git pull` there to get back to a clean tracked state.

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
