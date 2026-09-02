# CompareM2 v3 — decision log

How the design in [DESIGN.md](DESIGN.md) arrived. Read that first; this file is
for when you want to know *why* something is the way it is, or whether an idea
has already been tried and rejected.

**Rules for this file.** Entries are dated and append-only. If a decision is
reversed, edit the entry to say so and why — do not delete it, because the
reasoning that reversed it is usually worth more than the decision. Reversed
entries are the most useful thing here.

Measurements are labelled measured or documented. Never overwrite a measurement
with a guess.

---

## 2026-09-01

### Lives on branch `v3` in `cmkobel/comparem2`
Becomes CompareM2 v3 rather than a separate project, so the name, citation and
docs domain carry over. v2 was left in place on the branch for reference.

**Completed 2026-09-02** (`c7e8d91`): v2 removed — 92 files, `workflow/`,
`dynamic_report/`, `profile/`, `resources/`, `config/`, the launcher, Dockerfile,
environment.yaml, changelog.txt, the Rproj. History keeps all of it and `master`
still serves it. Deleting it forced rewrites of `pixi.toml`, CI, `CLAUDE.md`,
`README.md` and `docs/`, which had all been describing v2.

Repository history is **not** rewritten. Purging the 35.4 MB of PDFs
accidentally committed in `952a7d0` would only be worthwhile alongside the
157 MB of raw genome FASTA, and that sits below the branch point — so it would
rewrite `master` and every published SHA of v2 across ten remote branches. Not
worth the disk.

### Python only; R is dropped
CLI, TUI and report share one runtime. Removes R, pandoc, rmarkdown, tidyverse
and r-clusterProfiler — the heaviest non-database dependency in v2. The sixteen
existing `.rmd` report sections are not ported; they are rewritten.

### Snakemake is kept, as an executor rather than a framework
**Reversed within the day.** The first call was to drop Snakemake on the grounds
that a single-digit tool count makes a DAG engine overkill. That reasoning was
wrong: the DAG was never the hard part. The expected deployment is a workstation
or an HPC head node submitting to SLURM, and *cluster execution* is the hard
part — sbatch generation, `afterok` chains, `squeue`/`sacct` polling,
partial-failure recovery, retries. Snakemake's SLURM executor plugin already
solves it. The install-weight argument also fails on the numbers: Snakemake is
tens of MB against 141 GB of databases.

Cost accepted knowingly: **the TUI gets harder.** Snakemake owns the event loop
and its logging, so a live interface has to drive it through `snakemake.api` with
a custom log handler and depends on that log event schema. Resolved later the
same day — see below.

### Wide, not deep
See DESIGN.md. Consequence: antismash, gapseq, interproscan and eggnog are out
of the default path, which is what makes one environment possible.

### One environment, plus CheckM2 on its own — revised same day
**The first version of this decision was wrong, and how it was wrong is the
point.** The original solve succeeded, so "one environment" was recorded as
settled. It only succeeded because the solver quietly selected **bakta 1.8.1**
(2023), the newest Bakta that could co-exist with CheckM2 — and 1.8.1 then
crashes against pyrodigal 3.x, which renamed `OrfFinder` to `GeneFinder`.
**A solve is not a working environment**, and only running the tool revealed it.

Isolated afterwards, on `linux-64`:

| Combination | Result |
| ----------- | ------ |
| bakta ≥1.10 alone | solves — 1.12.1, diamond 2.2.5 |
| bakta ≥1.10 + checkm2 | **conflict** |
| bakta ≥1.10 + the other tools | solves |

CheckM2 pins DIAMOND 2.1.x; current Bakta needs 2.2.x. CheckM2 is too valuable
to drop, so it is the one tool marked `isolated=True`.

### The superseded single-environment measurement
Kept because it is the measurement the wrong decision above was based on.
Measured with pixi against `linux-64`:

| Environment | Packages | Download |
| ----------- | -------: | -------: |
| core 8 tools | 576 | 0.76 GB |
| core + bakta + snakemake | 830 | 0.83 GB |
| core + prokka + snakemake | 840 | 1.60 GB |

Prokka costs roughly twice bakta in package weight — it pulls perl and BLAST,
which is why bakta is the only annotator.

### The tool set: 14 in, 8 out
Selected tool by tool. v2 had 20. **Now thirteen — sylph removed 2026-09-02**,
see the post-mortem below.

**Kept** — seqkit, checkm2, gtdbtk, bakta, amrfinder, mlst, mashtree,
treecluster, skani, panaroo, snp-dists, fasttree, carveme.

**Dropped** — assembly-stats (seqkit covers it), prokka, eggnog, dbcan,
interproscan, gapseq, iqtree, clusterProfiler (an R package, and R is gone).

**New** — skani, for all-against-all ANI. sylph was also selected here and later
removed.

Two worth their reasoning:

- **carveme in, gapseq out.** Genome-scale metabolic models are the one
  capability Bactopia, funcscan and DRAM2 all lack. CarveMe keeps it at minutes
  rather than gapseq's hours, and solves with open-source SCIP — no CPLEX
  licence.
- **antismash dropped.** Selected, then dropped when it turned out to break the
  single environment: it pins `biopython 1.78` and `diamond 2.1.11` against the
  newer versions checkm2, bakta and gtdbtk require. One environment was judged
  worth more than one BGC caller, especially as funcscan ships four.

### Bakta light, not full
v2 downloads `--type full` (30 GB compressed, 84 GB on disk). v3 uses `light`.
Saves 29 GB for less specific annotation, which a wide view can absorb — but
note the Bakta paper's annotation-quality figures are measured on the full
database and do not transfer.

### Vertical slice runs end to end
seqkit → mashtree → treecluster on four *E. faecium* genomes, 7/7 steps, report
rendered. This closed the central architectural question: **Snakemake rules can
be generated from declarative specs**, including wildcards, dependencies and
stdout redirection. Three bugs it caught are in the post-mortems below.

### Passthrough parameters survive from v2
v2's `set_<tool>--<flag>: <value>` becomes `--set tool--flag=value`, backed by a
`params` field on each `Tool`. Defaults carried over from v2's `config.yaml`
verbatim, so v3 reproduces v2's behaviour unless told otherwise.

### Unit tests, which v2 never had
v2 validated only by running the pipeline in CI — the wrong instrument for a
generator. Started at 26 tests.

### TUI against Snakemake — no fork needed
Snakemake 9.26.1 ships a logger plugin system (`LogHandlerBase`, `LogEvent`,
`LoggerPluginRegistry`), and events can also be captured in-process by attaching
a `logging.Handler` to the `snakemake` logger while driving `SnakemakeApi`. A
spike captured 19 events on a two-job run, with structured fields rather than
scraped text: `run_info`, `job_info`, `job_started`/`job_finished`, `progress`,
`shellcmd`, `job_error`/`error`.

That is everything a live keyboard-driven interface needs. The offer to ship a
modified Snakemake stays unused, which is the better outcome — a fork would have
to be rebased forever.

---

## 2026-09-02

### Pin a minimum version for every tool
See the post-mortem below. The rule: pin a minimum in each spec's `conda` field,
and treat verification as tracking *execution*, never installation.

### The report explains itself, from the tools' own papers
The CompareM2 paper claims the report is accessible to non-bioinformaticians.
v2 delivered that with hand-written RMarkdown per section; v3 had one-line
`summary` strings and nothing else.

Every tool now carries a `Guidance` value — see DESIGN.md for the shape. Three
decisions inside it:

- **Guidance lives outside `catalogue.py`**, with a test enforcing completeness
  rather than proximity.
- **Collapsed by default.** An expert should not scroll past a page of prose to
  reach their data; a non-expert needs it one click away rather than in a manual.
- **Every number is quoted and grep-checked.** 178 quantitative claims were
  extracted with a verbatim quote each and verified by substring match against
  `pdftotext` output; all 178 passed. Where a statement is ordinary
  methodological caution rather than a paper's finding, the sentence says so.

The 24 papers are in `papers/` (untracked) with `papers/tools.bib`, and
`papers/SUMMARIES.md` holds the long-form reading including what the papers do
*not* answer.

### The unit tests are in the repository
`.gitignore` carried a blanket `tests` rule whose only live effect was keeping
`tests/unit/test_v3.py` out of git — the suite presented as v3's main improvement
over v2 existed in one working tree and could not run in CI. Nothing else under
`tests/` was caught by it. Replaced with narrow rules for unpacked genomes and
run output.

### Two docs pages are generated from the specs
`docs/30 what analyses does it do.md` and `docs/99 citation.md` are written by
`docs/generate.py`, and CI fails if they are stale. v2 hand-maintained both and
both drifted: the citation list kept dropped tools and missed added ones, and the
analyses page documented defaults that no longer matched `config.yaml`.

### Databases download themselves, as Snakemake rules
See DESIGN.md for the mechanism and the post-mortem below for what was wrong.
The fetch mechanisms are not uniform, which is presumably why it was deferred:

| Database | Fetch | Ready |
| -------- | ----- | ----- |
| checkm2 | curl + tar, then link the release's `.dmnd` to a stable name | `checkm2/checkm2.dmnd` |
| bakta-light | `bakta_db download`, then move its `db-light` into place | `bakta/version.json` |
| gtdb | curl + tar `--strip-components=1` | `gtdb/.fetched` (stamp) |
| amrfinder | `amrfinder -u` | `amrfinder/.updated` (stamp) |

Two use stamps and both reasons are measured, not assumed. **amrfinder**:
`amrfinder -u -d <dir>` exits with *"AMRFinder update option (-u/--update) only
operates on the default database directory. The -d/--database option is not
permitted"* — so its data goes into `$CONDA_PREFIX` and cannot live under
`--databases`. `Database.out_of_tree` records that rather than letting the spec
imply otherwise. **gtdb**: 141.4 GB has never been downloaded, so no interior
filename can be asserted; a stamp is the alternative to inventing one.

### Working checkout moved off scratch
Verification ran for a day in `/evo/postdoc/cm2v3`, a gitignored scratch
directory inside the `postdoc` repo, synced by rsync — which would silently
drift from the branch. Replaced by a real clone. Current paths are in
[STATUS.md](STATUS.md).

Two operational facts worth keeping: moving a pixi project invalidates its
environments, because conda bakes the absolute prefix into shebangs and RPATHs,
so `rm -rf .pixi && pixi install` is required after any move. And GenomeDK was
the alternative machine and is not reachable non-interactively
(`Permission denied (publickey,keyboard-interactive)`).

---

## What went wrong

Each of these produced a rule in DESIGN.md. They are collected here so the
decisions above stay readable, but the rule and its reason belong together —
if you are about to undo one of those rules, this is the section to read first.

### A solve is not a working environment
Two tools silently resolved to years-old builds that install cleanly and crash on
first use:

| Tool | Unpinned resolved to | Breaks on |
| ---- | -------------------- | --------- |
| bakta | 1.8.1 (2023) | `pyrodigal.OrfFinder`, renamed `GeneFinder` in pyrodigal 3.x |
| panaroo | 1.1.2 (2020) | `Bio.Alphabet`, removed in Biopython 1.78 |

Neither is a solver failure — the solver did what it was asked. Bioconda keeps
old builds forever, and an unconstrained `tool = "*"` lets the solver satisfy
some *other* package's constraint by reaching back years.

### Three bugs the first end-to-end slice caught
None of which reading could have found:

1. **Bare `{sample}` in a shell block is a runtime NameError.** Snakemake
   resolves wildcards in `input:`/`output:` but exposes them as
   `wildcards.sample` inside `shell:`. The Snakefile parsed and built the right
   DAG, then failed on execution.
2. **`TreeCluster.py` requires `-t/--threshold`.** The drafted command omitted
   it. v2's defaults (`max_clade`, `0.05`) were adopted.
3. **Sample names can contain spaces.** v2's own test set ships
   `116_2 duplicate.fna`, which produced a silently broken wildcard.

### Four defects found by reading the tools' papers
Writing interpretation guidance meant checking each declared command against
what its paper says the tool does. Every one of these produces a Snakefile that
parses and a DAG that builds, which is why no test caught them:

| Tool | Defect | Fix |
| ---- | ------ | --- |
| sylph | `sylph profile <db>.syldb <assembly>.fna` passes assemblies as positional FASTA, which `profile` treats as **reference genomes**, not samples. With no FASTQ or `.sylsp` present it exits `No read files found`. | **removed** |
| skani | Ran at defaults `k=15, c=125`, described in the paper as tuned for complete, similar genomes. | `-c 70`, the paper's middle preset |
| panaroo | The report's Soft core (95–99%) and Cloud (<15%) bins are unreachable below 20 and 7 genomes, so two rows always read 0 and all accessory content piled into Shell. | exact counts below 20 genomes |
| fasttree | Declared `threads=4` but invoked plain `FastTree`, which is single-threaded. | `threads=1` |

**sylph was removed, not repaired, and that is the lesson.** It profiles
metagenomic *reads*; v3's input is assemblies, which is not the question it
answers. The broken command line was a symptom of selecting a tool by capability
blurb rather than by input type. skani already covers fast
assembly-to-assembly identity, so nothing was lost.

**The panaroo bin arithmetic was real on live data, not just in theory.** On four
genomes the old bins would have read `Core 2091 / Soft core 0 / Shell 1689 /
Cloud 0`. The switch-over is now pinned by tests, including one asserting that 20
is the smallest N at which all four conventional bins can hold a cluster.

### `Database.url` was dead code
`Database` carried a `url` for all four databases and **no code read it** — the
only `.url` in `src/` was `Citation.url`. So v3 declared four databases with
measured sizes, printed `databases: 143.2 GB + 2 of unknown size`, and then
never fetched anything. That is the worst of the three options: not having
downloads is defensible, announcing a total and failing minutes later inside a
tool is not.

It went unnoticed because the test machine's databases had been placed by hand.

### GTDB-Tk's database was unreachable through `--databases`
GTDB-Tk has no flag for its database; it reads `GTDBTK_DATA_PATH`, which the
generated rule set nowhere. So `--databases` was **silently ignored for the
largest database in the pipeline** — 91% of the install weight, pointed wherever
the ambient environment happened to point. `Tool.env` exists because of this.

### Snakemake locks the working directory, not the output directory
Two `cm2` runs in one checkout collided even with different `--output`, and a
killed run left a lock the next could not clear. The cause: Snakemake locks its
*working directory*, and `cli.py` ran it from the checkout root while pointing
`--snakefile` into the output directory.

Fixed by resolving `--output` to an absolute path and passing it as
`--directory`, so `.snakemake` lives in the output directory and the lock means
what it should. Rejected `--nolock`, which removes the guard entirely and would
let two concurrent runs on one output silently clobber each other. Cost,
accepted: absolute paths in the generated Snakefile, which was already true of
the input symlinks.

### `--set` could not express a single-dash flag
Found by adding skani's `-c 70`, the first single-dash parameter in the
catalogue. `parse_overrides` split on the first `--` and prepended `--` to the
rest:

| Written | Produced | |
| ------- | -------- | - |
| `skani-c=125` | `SystemExit` | rejected outright |
| `skani--c=125` | `-c 70 --c 125` | default not replaced, and `--c` is not a skani flag |

The flag now keeps whatever dashes it was written with, and the tool is matched
against the catalogue longest-first so `snp-dists` is not split on its own
hyphen. Every other tool's parameters happen to use long flags, which is why it
survived — a passthrough mechanism only gets exercised where a default exists to
exercise it.

### skani was throwing away half its output
`skani triangle --full-matrix` writes two matrices: identity at `-o`, and the
aligned fraction at that path plus `.af`. v3 declared only the first. Found by
reading a real run's log, not the paper — the paper says AF is computed, not
that the file appears next to the output.

It matters because skani emits an ANI once alignment covers as little as ~15% of
a genome, so identity alone cannot distinguish whole-genome relatedness from a
shared plasmid or conserved core. The four test genomes show it: 116_2 against
E8202 reads 99.14% ANI on an aligned fraction of 74–90%, depending on direction.

### A comment that was simply wrong
`pixi.toml` claimed the test genomes are listed rather than globbed because a
shell glob would split `116_2 duplicate.fna` into two arguments. False — glob
results are not word-split, only unquoted variable expansion is. Corrected after
a real run expanded `tests/E._faecium/*.fna` to four correct paths.

Kept as an entry because it is the same failure mode as the rest of this
section, applied to prose instead of code: a plausible mechanism asserted
without being tested.
