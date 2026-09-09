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

**Amended 2026-09-02.** The spike counted events and did not read their fields.
`job_finished` turned out to carry neither `rule_name` nor `jobid` — it carries
`job_id`, with an underscore — so "the events exist" was true and "the events
say which rule finished" was not. See *The first `--tui` run* below.

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

### v3 becomes `master`, and v2 is not preserved on a branch
Merged as a **fast-forward** — `v3` was 24 commits ahead of `master` and 0
behind, so no merge commit and no conflict. The default branch on GitHub was
already `master`, so nothing had to be repointed.

**v2 is deliberately not given a branch.** Carl's call, and worth recording
because the pre-merge survey argued the other way. Three facts that made the
argument, all measured:

- `origin/v2` is a 2019-era branch, **2,008 commits** behind the pre-merge
  `master`. It is *not* what the paper describes, despite an earlier claim in
  STATUS.md that said so — that claim was wrong and is corrected.
- The last v2 tag, `v2.9.1`, is **443 commits** behind the pre-merge `master`
  tip. Those 443 include the Bioinformatics paper, the Snakemake 7→9 migration
  and the docs rewrite, so no existing tag names the state v2 actually ended in.
- That tip survives in history and on the `ai-1` branch. Reachable, just not
  under an obvious name.

The consequence to accept: reconstructing a runnable v2 means finding a SHA, not
checking out a branch. The v2 *documentation* is unaffected — Read the Docs
builds `stable` from the newest tag, which is still `v2.9.1`, so
`/en/stable/` keeps serving v2 while `/en/latest/` follows `master` to v3.

### Docs stop describing v3 as a side branch
The merge invalidated the framing rather than any code. README and
`docs/index.md` told the reader to use `master` for the real thing — which is
now v3 itself — and `docs/05` and `docs/10` both said `git clone -b v3`. All
reframed to *pre-release*: v3 is the default branch, and what is missing is a
bioconda package and a container image, not a branch to find it on.

Version stays `3.0.0.dev0`. **Merging is not releasing**, and the release story
still needs the bioconda recipe, a container, and GTDB-Tk actually executed.

Also corrected here, found while editing: the unit-test count was given as 75 in
two places and 100 in a third, against 105 collected (103 passing, 2 skipped
without Snakemake, 1.05 s). And `docs/10` pointed at `DESIGN.md` for the
execution-status table, which lives in `STATUS.md`.

### The docs requirements are floors, not a lock
Publishing v3 to `master` moved four open Dependabot alerts onto the default
branch — 1 high, 2 moderate, 1 low, all in `docs/requirements.txt`, which was a
full transitive lock compiled in 2022 against Python 3.10 and never revisited.

Three of the four came from `pymdown-extensions`, pulled in only by
`mkdocstrings[python]`, which nothing used: `mkdocs.yml` loads the `search`
plugin alone and no page has a `:::` autodoc directive. Removed rather than
bumped. The fourth was `Markdown` 3.3.7, which floors resolve to 3.9.

Replaced the lock with `mkdocs>=1.6` and `markdown-include>=0.8`. A docs build
is not what a lock buys much for — `mkdocs build --strict` fails loudly and
immediately — and the lock's real effect was recurring alert noise. Verified in
a clean virtualenv: those two floors plus Markdown 3.9, no pymdown-extensions,
`mkdocs build --strict` clean.

### Relative paths resolve against `$INIT_CWD`, not the cwd
`pixi run cm2 *.fna --tui`, typed in `tests/E._faecium/`, reported
`no such file: 116_2 duplicate.fna, 116_2.fna, E8202.fna, SRR24.fna` — all four
of them sitting in the `ls` directly above it. Not the space in the filename,
which survives intact: **a pixi task executes from the workspace manifest root,
not from the shell's directory**, so the four relative names were looked up one
directory tree away.

Measured against pixi 0.78.0, from `<root>/sub` containing `a b.txt` and
`c.txt`:

```
✨ Pixi task (count): python3 -c '...' a b.txt c.txt
2 ['a b.txt', 'c.txt']      # two arguments, the space intact
PWD=/private/tmp/pixicwd            # the manifest root
INIT_CWD=/private/tmp/pixicwd/sub   # where it was typed
```

So `cli.py` now resolves every user-supplied path — inputs, `--output`,
`--databases` — against `$INIT_CWD` when it names a real directory, and the cwd
otherwise. Absolute paths are untouched; `base / path` is a no-op on them.

Rejected: a `cwd` on the task, which is static and cannot mean "wherever the
user is"; and telling users to type `"$PWD"/*.fna`, which is a workaround
written in the documentation rather than a fix, and leaves `--output` still
defaulting to a `results_comparem2` beside the manifest instead of beside the
genomes.

The error message now names the resolved path. `no such file: 116_2.fna` is
least useful in exactly the case that produced it — when the file is in the
directory the user is looking at — and the absolute form says where the lookup
happened.

### The remote follows the repository rename
`origin` now points at `cmkobel/CompareM2`, and the `cmkobel/comparem2` spelling
is gone from `mkdocs.yml` and the docs. GitHub redirected the old name, so this
fixes a notice on every push rather than a breakage. Occurrences inside dated
entries in this file are left as written — the log records what was true then.

### The bioconda package ships the pipeline, not the tools
The choice was between two recipes, and Carl took the first:

**A. Pipeline only.** Run dependencies are `python`, `snakemake-minimal`,
`textual`, the two executor plugins, and `conda`. Every analysis tool arrives
through `--software-deployment-method conda`, deployed from the `envs/*.yaml`
that `prepare()` already writes.

**B. The tools as run dependencies.** Rejected. It cannot include all thirteen —
CheckM2 pins DIAMOND 2.1.x against Bakta's 2.2.x — so it would ship twelve and
handle the thirteenth some other way, put a thirteen-tool solve inside
bioconda's CI, and break whenever any of the thirteen changed upstream.

This is the model v2 used too (25 environments deployed at run time, against
v3's 14), and it is why the recipe for a thirteen-tool pipeline is 60 lines.

What it cost in code:

- `pyproject.toml`, which did not exist. `pixi.toml` said outright that this
  was "a workflow application, not a distributable package" — that sentence was
  the blocker. setuptools, `dynamic` version read from `__init__.py`, and two
  entry points (`comparem2`, `cm2`).
- `--use-conda` and `--conda-prefix`, wired through both execution paths: the
  CLI's Snakemake subprocess and the TUI's `SnakemakeApi` call, which needed
  `DeploymentSettings(deployment_method={DeploymentMethod.CONDA})`.
- `Database.conda`. A download rule is a rule and needs an environment like any
  other, and two of the four fetches run a tool binary rather than curl.
- A PATH preflight. A fresh `conda install comparem2` has none of the thirteen
  tools, and without this the first thing a new user saw was a Snakemake
  traceback from whichever rule was scheduled first.

Version stays `3.0.0.dev0`; `recipe/README.md` holds the release steps.

### Two fields for the steps around a command
GTDB-Tk needs a file written before it runs and its output reshaped after, and
neither is an argument list. Rather than let one tool become a hand-written
rule, `Tool` gained two fields:

- **`files`** — path-to-content, rendered from the `Context`, written by
  `prepare()` and declared as a rule input. GTDB-Tk's two-column `--batchfile`,
  which exists because canonicalisation puts each genome in its own directory
  and `--genome_dir` therefore cannot be used.
- **`post`** — argument lists run after the command. The `bac120`/`ar53` merge.

Considered and rejected:

- **`awk`/`printf` in the rule.** Untestable, and the quoting of a tab-separated
  format inside a generated shell block is exactly the kind of thing that looks
  right and is not.
- **A pre-step symmetrical to `post`.** Unnecessary: the batchfile's content
  depends only on the sample list, which is known before anything runs, so it is
  data rather than a step. Fewer moving parts, and it becomes a declared input,
  which a shell step never could.
- **Declaring GTDB-Tk's real outputs instead of merging.** An all-bacterial set
  writes no `ar53` file, so a declared output would fail on most real input.

Two details that took a moment to see. The step runs through an absolute
`sys.executable`, because under `--use-conda` the rule's environment holds the
tool and not CompareM2. And a declared file is written **only when its content
changes** — it is a rule input, so rewriting an identical one moves its mtime
and re-runs GTDB-Tk, hours of work triggered by a four-line file.

A test asserts GTDB-Tk is the only tool using either field. `post` is not a
licence for shell work.

### No container image
Dropped as a release blocker, on Carl's call: pixi and conda have both got good
enough that a hand-built image would be a third installation path to keep in
step with the other two, for a case neither of them already covers badly.

Bioconda builds a BioContainer for the package automatically, so an image will
exist regardless — it just contains the pipeline and no analysis tools, which
is honest about what the package is. A thirteen-tool image was never going to
be one environment anyway (CheckM2's DIAMOND 2.1.x against Bakta's 2.2.x), so
it would have been two environments in one image, built and verified by hand.

v2 had one. This reverses "both existed for v2 and will return", which was
written in the README and in STATUS.md as recently as this morning.

### Environments are addressed by content, which AMRFinder needs
Read from Snakemake 9.26.1's `deployment/conda.py`: a deployed environment's
directory is `md5(realpath(envs_dir) + env file content)`, deliberately
excluding the env file's own path. Two consequences that are not cosmetic.

**It is what makes AMRFinder work at all under `--use-conda`.** `amrfinder -u`
refuses `-d` and writes into `$CONDA_PREFIX`, so the download rule and the
analysis rules have to end up in the *same* deployed environment — which means
their env files must be byte-identical, not merely equivalent. Hence
`_AMRFINDER_SPEC` as a shared constant, and a test asserting the two rendered
files are equal.

**It is why `--conda-prefix` defaults to a shared `~/.comparem2/envs`.** Moving
the prefix changes every hash, so all 14 environments re-solve and AMRFinder's
database is fetched again. Same reasoning as the database directory, one step
sharper.

### Snakemake by name, not by module
Found by running the installed package: `comparem2 ... --use-conda --dry-run`
from `/tmp/cm2venv/bin/comparem2`, without the environment activated, died with
`FileNotFoundError: 'snakemake'` out of `subprocess.run` — three lines after
printing what it was about to do. The CLI shelled out to a bare `snakemake`,
which is only on PATH when the environment happens to be active. It is now
`sys.executable -m snakemake`: this package's own dependency, so the right copy
is the one beside the running interpreter.

Under pixi this never surfaced, because `pixi run` activates the environment.

### CarveMe gets a wrapper, because a solver parameter is not reachable
`carve` spent 601 s of a 609 s run inside SCIP and then returned a model with
253 of its annotated reactions missing, because conda-forge's SCIP links PaPILO
and PaPILO presolves this problem's optimum away (see *What went wrong*). One
SCIP parameter fixes it. Neither CarveMe nor ReFramed exposes one, so
`src/comparem2/carve_scip.py` sets it and calls CarveMe in-process.

Considered and rejected:

- **Relaxing CarveMe's gap limit.** The cheapest-looking knob, and wrong: at a
  5% gap the run finishes in 33 s with a *1,160-reaction* model, worse than
  either. It buys speed by keeping the bad search's first feasible point.
  Measured, not reasoned about.
- **Pinning `scip`.** conda-forge has 10.0.2, and it carries PaPILO 3.0.0 and
  behaves the same — still in the MILP at 300 s. There is no conda-forge SCIP
  without PaPILO to pin to.
- **`pyscipopt` from PyPI**, whose wheel bundles a SCIP without PaPILO and
  solves in 9.8 s. It would mean a pip package shadowing the conda one that
  bioconda's carveme depends on, in an environment pixi and conda both manage —
  and nothing like it would work under `--use-conda`.
- **A rewrite, in Rust or otherwise.** The question that started this. Every
  line of Python in the path — CarveMe's included — is about 8 s of the 609:
  DIAMOND 1.8 s, universe model 1.7 s, scoring 4.0 s, SBML out 0.9 s. A perfect
  rewrite of all of it saves under 2%.
- **More threads.** SCIP's branch-and-bound here is one core, and CarveMe is
  per-genome, so the fan-out that matters is Snakemake's and already exists.

The wrapper runs under a bare `python`, which contradicts the rule `steps.py`
established a few entries above — and has to. `steps.py` is our code and needs
our interpreter; this imports `carveme` and needs the interpreter of the
environment the tool is in, which under `--use-conda` is a different one. Both
import nothing from their own package, for the same reason, and a test enforces
both halves.

It also earns its keep twice: giving `carve` an input path inside CarveMe's own
output directory is what stops it overwriting Bakta's feature table.
`Tool.executable` came with it, because the preflight now sees an interpreter as
argv[0] and would have looked for `python`, found it, and never reported carveme
missing.

---

## 2026-09-03

### FastTreeMP was reachable, measured, and rejected
The spec's comment claimed `threads=1` because there was "nowhere to set"
`OMP_NUM_THREADS` — commands are argument lists, not shell strings. That reason
had already expired: `Tool.env` was added for GTDB-Tk's `GTDBTK_DATA_PATH` and
`snakefile.py` emits an `export` line for it, and `bioconda::fasttree` ships
`FastTreeMP` in the same package as `FastTree`. The switch was a two-line change.

**Measured** on thylakoid, `fasttree 2.2.0` build `h7b50bb2_1` — the build in
`pixi.lock` — with the catalogue's own flags, `-nt -gtr`.

A real panaroo core alignment, 7 taxa x 2,066,459 bp:

| binary | threads | wall s |
| --- | ---: | ---: |
| FastTree | 1 | 242.3, 209.4 (two reps) |
| FastTreeMP | 1 | 251.7 |
| FastTreeMP | 2 | 278.5 |
| FastTreeMP | 4 | 254.6, 200.2 |
| FastTreeMP | 8 | 204.3 |
| FastTreeMP | 16 | 208.2 |

Two reps of the *same single-threaded binary* differ by 33 s, 14% — wider than
any plain-against-MP gap in the table, and the machine was carrying load 5–9
from other work throughout. On the quiet 4-taxon case (`tests/E._faecium`,
1,935,176 bp) the variance is small and the sign is consistent: 22.67 s plain
against 23.19 s at 4 threads and 23.23 s at 8, so MP is 2.3–2.5% **slower**.

Simulated alignments were needed to vary the one axis real data here cannot —
every core alignment on the machine has 4 or 7 taxa. 100 kbp, ~0.5% divergence
per branch:

| taxa | plain | MP x8 | speedup | CPU s, plain → MP |
| ---: | ---: | ---: | ---: | --- |
| 25 | 25.13 | 22.31 | 1.13x | 25.1 → 35.9 |
| 100 | 132.54 | 115.67 | 1.15x | 132.5 → 181.9 |
| 400 | 597.47 | 546.18 | 1.09x | 597.5 → 854.9 |

**The speedup does not grow with taxon count.** That was the obvious objection —
four taxa give FastTree one internal split to spread across threads, so of
course OpenMP does nothing; hundreds of genomes should be different. They are
not. From 25 to 400 taxa the payoff sits between 1.09x and 1.15x, for 40–65%
more CPU and 8 cores Snakemake would have reserved.

The 25- and 100-taxon runs are single measurements on a quiet machine (load
1.2–2.0). The 400-taxon row is the mean of two runs each, ordered
plain/MP/MP/plain so that drift in machine load cancels: the machine was busy
with other sessions and load fell monotonically across the four runs, 10.29 →
9.70 → 8.63 → 6.28 sampled during each. The plain runs therefore averaged load
8.29 against the MP runs' 9.17, which biases *against* MP — 1.09x is a floor,
not a ceiling.

**A first pass at 400 taxa was thrown away, and it is worth saying why.** It
reported MP at 580.9 s and 589.1 s against a plain baseline of 545.6 s, i.e. MP
slower, and that was an artifact: load was sampled once before each run rather
than throughout, and the baseline happened to start at load 1.10 while the MP
runs started at 7.11 and 9.99. Re-running the plain binary under realistic load
gave 615.1 s for the same work — a 13% penalty that had been credited to
FastTreeMP. The conclusion did not change, but the number did, and an
unbalanced benchmark on a shared machine is how you get a right answer for a
wrong reason.

Wall time equalled user+sys in every plain run. Every MP run exceeded it by
1.5–2.0x, whether given 2 threads or 16 — that is the average number of cores
actually kept busy, and it is what caps the speedup. It is not Amdahl on the
support phase either: with `-nosupport` removing that stage entirely, MP at 16
threads was still slower than one thread (619.0 s against 526.2 s).

One side effect, had it been adopted: at 7 taxa MP changes branch lengths in the
sixth decimal — `0.003849049` against `0.003847998` — from floating-point
summation order. Topology and support values are identical, MP is deterministic
across thread counts, and at 4 taxa the output is byte-identical.

### 59% of FastTree's runtime is support values nothing reads
Found while reading the phase log for the above. At 7 taxa the run reaches
"Optimize all lengths" at 99.25 s and reports `Total time: 239.96` — everything
after is SH-like support, 1000 resamples. Confirmed by running the flag:

| input | default | `-nosupport` | saved |
| --- | ---: | ---: | ---: |
| 7 taxa x 2.07 Mbp | 239.96 s | 98.03 s | 141.9 s, 59% |
| 4 taxa x 1.94 Mbp | 21.46 s | 12.93 s | 8.5 s, 40% |
| 400 taxa x 100 kbp | 545.58 s | 526.24 s | 19.3 s, 3.5% |

The fraction is a property of the shape, not of FastTree: support cost scales
with alignment length, the ML search with taxon count. Few genomes over a
megabase core alignment is the pipeline's normal case, so 59% is the figure
that applies here.

Topology and branch lengths come out byte-identical; only the `1.000` internal
labels drop. And `report.py`'s `draw_tree` writes text in an `elif node.name`
hanging off `if node.children`, so it labels leaves only — FastTree's support
values are parsed into `_Node.name` and never drawn. At the 4-taxon test-set
shape they do not even reach the newick: with 3 unique sequences of 4 there is
no non-trivial split to label, and `n4_plain.newick` and `n4_nosup.newick` are
byte-identical, so those 8.5 s produce nothing at all.

Not acted on, because `fasttree.newick` is a product in its own right and
someone opening it in iTOL at larger taxon counts would lose real information.
The choice is to drop the values or to render them, not to leave both.

### A fourteenth tool, and the niche-media idea that did not survive measurement
The ask was to make `carveme` more useful: simulate growth on media standing in
for ecological niches, score each genome on them, and get a high-level view of
its metabolism. The idea is right and the implementation it implies does not
work. What shipped instead is `biosynthesis`, per-compound rather than
per-medium.

**No draft grows on any defined medium.** Eleven real CarveMe models — four
*E. faecium* from `results_full13`, seven *S. aureus* from the `verify` run —
against CarveMe's own media, applied exclusively (every exchange to 0, then the
medium's compounds to −10 mmol/gDW/h, which is the paper's phenotype-array
protocol):

| medium | growing, of 11 |
| --- | ---: |
| M9 | 0 |
| M9[-O2] | 0 |
| M9[glyc] | 0 |
| LB | 0 |
| LB[-O2] | 0 |
| complete, every exchange at −10 | 11, at 9.9–21.1 h⁻¹ |

Cross-checked in two libraries — cobrapy 0.32.1 with GLPK and reframed 1.6.0
with SCIP — which disagree on the complete-medium *rate* and agree on every
zero. The rates differ because cobrapy's SBML reader adds exchange reactions for
`boundaryCondition="true"` species, of which CarveMe marks every extracellular
metabolite, so it solves a model with six more exchanges than CarveMe wrote.
ReFramed reads what is there, which is why the tool uses it.

**Growth is one bit that a single metabolite destroys.** Adding a demand
reaction per biomass precursor: on LB, `116_2` can produce 52 of its 53
precursors and fails on menaquinol-8. `COL` fails on asparagine, alone. Per
compound, 52 bits survive what kills the one — which is the whole argument for
the shape the tool ended up with.

**A plain M9 producibility scan cascades**, so it cannot be the readout either:
`116_2` comes out with 26 of 53 precursors unreachable, but no folate means no
purines means no ATP, and reading 26 as 26 auxotrophies is wrong. Hence three
verdicts — `de_novo` from M9 alone, `upstream` from M9 plus the rest of the
panel, `none` from neither.

**The background may not be the complete medium.** With every exchange open,
`COL` looks able to make asparagine: it imports the Gly-Asn dipeptide and
hydrolyses it. Salvage is not synthesis, and the panel background excludes
peptides by construction because they are not on the panel.

Two other framings were built and rejected on the way:

- *Greedy medium reduction* — close each non-M9 exchange in turn, keep it closed
  while growth survives. Fast (211 LPs, 1.4 s) and it separates the species
  cleanly, but it is order-dependent and says so out loud: it reports `116_2` as
  requiring L-methionine *S*-oxide (because `met__L` sorts first and was closed
  first), the tripeptide Pro-His-Glu, and benzoate. Nobody can act on that.
- *A Biolog-style carbon-source scan*, the one framing with published
  validation. It works — 196 substrates, 3 s a genome, and it recovers real
  biology: the *E. faecium* drafts use sugars and no TCA intermediates, which is
  correct for an organism without a complete TCA cycle, and the *S. aureus*
  drafts use αKG, succinate, fumarate and malate. But it is not comparable
  across genomes. It needs a background found by search, and `E8202` returned 11
  usable substrates against 53 and 67 for its two conspecifics purely because
  its background leaks 0.0154 h⁻¹ against their 0.0064 — the score tracks the
  background, not the metabolism. The panel size varies too, by 40–60 substrates
  between strains of one species, because a model only carries exchanges for
  what carving kept.

### The biosynthesis probe target is the pathway's product, not the vitamin
Three targets were tried and rejected, each against `iML1515` — *E. coli* K-12,
manually curated — where the answer is known:

| probed | verdict on iML1515 | why it is the wrong target |
| --- | --- | --- |
| `fol` | no route | de novo synthesis runs dihydropteroate → dihydrofolate → THF; folate is not an intermediate. It also called *S. aureus* unable to make the compound sulfonamides work by blocking. Use `thf` |
| `thm` | no route | free thiamine is a salvage substrate; synthesis ends at thiamine phosphate → ThDP. It called *E. coli* a thiamine auxotroph. Use `thmpp` |
| `lipoate` | no route | the de novo product is protein-bound lipoyl, not free lipoate. Returned "cannot make" for all eleven drafts, and would for any organism. Dropped |

`nad` replaced `nac` for the same reason in reverse: with both on the panel each
rescued the other, so the pair measured a kinase. The rule that came out of it —
one member per nutrient family, probed in the form the pathway ends at — is in
DESIGN.md and has a test.

With those fixed, **iML1515 returns 31 of 32 de novo**. The exception is
adenosylcobalamin, which *E. coli* genuinely cannot synthesise de novo. That is
the calibration for the whole section.

**The family rule was then checked empirically rather than argued.** Leave one
compound out of the panel, recompute, and see whether any *other* compound's
verdict moves, over five models. 2 of 32 move, and both are real precursor
relationships rather than interconvertible forms:

- dropping `met__L` flips `ile__L` from `upstream` to `none` in `116_2` and
  `E8202` — isoleucine via 2-oxobutanoate from methionine, an alternative to the
  threonine route.
- dropping `thr__L` flips `gly` the same way in the same two — glycine from
  threonine by threonine aldolase, alongside the serine route.

Which is exactly what `upstream` is for. It also settled the one family the unit
test got wrong on first writing: protoheme and siroheme are separate branch
products of uroporphyrinogen III, not forms of one nutrient, and dropping either
changes no verdict in any of the five models.

### What the panel does and does not recover
Measured on the eleven drafts. `none` in all four *E. faecium*: leucine,
methionine, threonine, tryptophan, valine, riboflavin, pantothenate, NAD,
biotin, and both quinones — with arginine and histidine in three of the four,
`E8202` having them as `upstream`. `none` in all seven *S. aureus*: thiamine
diphosphate, NAD, biotin — and asparagine.

The B-vitamin requirements are right; enterococci were the assay organisms for
those vitamins historically. Thiamine and nicotinate are in the described
chemically defined medium for *S. aureus*. Menaquinone-8 comes out de novo in
four *S. aureus* and unreachable in every *E. faecium*, and ubiquinone-8 in
neither — correct for Firmicutes.

**Asparagine in all seven *S. aureus* is a false call**, and it is left in the
report rather than special-cased: a uniform column is flagged as carrying no
comparative information, which is the general form of the warning, and a panel
member removed because one clade gets it wrong stops being comparable.

The duplicate pair is identical on both tables, which is the standing
cross-check.

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

### GTDB-Tk's rule was wrong in four ways, and its comments said otherwise
The last unexecuted tool, taken apart before spending six hours downloading its
database. Every one of these would have surfaced only at runtime:

1. **The batchfile was never written.** `--batchfile` named a path no rule
   created. The comment above the command read "the rule writes a batchfile
   first".
2. **The summaries were never merged.** The declared `gtdbtk.summary.tsv` is not
   a name GTDB-Tk writes; it emits `bac120` and `ar53` separately. The comment
   read "the rule concatenates them".
3. **`--skip_ani_screen` does not exist in 2.7.2.** It was real in the versions
   whose ANI screen used Mash and needed either `--mash_db` or permission to
   skip. 2.7 screens with skani from the reference package and removed the flag,
   so the command would have exited on `unrecognized arguments`.
4. **The database was the wrong release.** The catalogue pointed at r226;
   `gtdbtk/config/common.py` in 2.7.2 reads
   `COMPATIBLE_REF_DATA_VERSIONS = ['r232']`. The tool would have refused the
   data after downloading 141.4 GB of it.

The pattern is the one this section keeps recording, in its purest form: two
comments describing steps that did not exist, in a rule nobody had run, plus
two facts about an installed tool that were never checked against the installed
tool. `gtdbtk --help` and one `grep` in the package answered both, in a minute.

Worth the note: (4) made things *better*. r232 is 60.8 GB against r226's
141.4 GB, because it dropped FastANI's reference genomes for skani sketches —
so the pipeline's largest cost more than halved, and the fix propagated to about
twenty numbers across the docs. `gtdbtk` is now pinned `>=2.7`, since the pin
and the database release are one decision, not two.

### A comment that was simply wrong
`pixi.toml` claimed the test genomes are listed rather than globbed because a
shell glob would split `116_2 duplicate.fna` into two arguments. False — glob
results are not word-split, only unquoted variable expansion is. Corrected after
a real run expanded `tests/E._faecium/*.fna` to four correct paths.

Kept as an entry because it is the same failure mode as the rest of this
section, applied to prose instead of code: a plausible mechanism asserted
without being tested.

### The first `--tui` run found five defects
Launched for the first time on thylakoid over the four *E. faecium* genomes.
checkm2 failed, and then the whole 13-tool run stopped: twelve tools that had
nothing to do with checkm2 stayed at `pending`, and a green `Report:` line was
printed for a run that had produced no output.

All five had one cause. `tui.execute()` duplicated `cli.py` instead of sharing
it, and drifted:

1. **The env files were never written.** `cli.py` called `render_envs()`; the
   TUI did not. checkm2's rule carries a `conda:` directive, so Snakemake ran
   the job and then failed *recording metadata* for it with a `WorkflowError`
   — which aborts the workflow, not the job. A missing 200-byte file killed a
   13-tool run.
2. **`--isolated-launcher` was dropped**, so `checkm2: command not found`.
3. **The Snakemake lock went back to the checkout**, because `runner.run()` was
   called without `workdir=`. Exactly the defect *Snakemake locks the working
   directory* records as fixed — fixed in the CLI, reintroduced beside it.
4. **`--until`, `--set`, `--keep-going` and `--dry-run` were ignored.** The TUI
   opened with all 13 tools selected regardless, putting GTDB-Tk's 141.4 GB one
   keypress away for a user who had asked for one tool.
5. **`runner.run(dry_run=True)` could not work.** `ExecutionSettings(dryrun=…)`
   raises `TypeError` — dry run is an executor plugin, not a setting — and the
   broad `except Exception` reported that as a failed workflow.

The fix is `snakefile.prepare()`, which writes the Snakefile *and* the env files
and is now the only way either entry point builds a workflow. The rest is
argument threading.

Two more surfaced while verifying the fix, both older than the run:

- **No tool could ever have been marked `done`.** `job_finished` carries
  `job_id`, and the runner read `jobid`, which is what `job_info` uses. Every
  finish event arrived with no rule attached and the TUI dropped it. A run that
  completed 6 of 6 steps reported that nothing had run. Only `job_info` carries
  `rule_name`, so the runner now remembers the mapping. Invisible before,
  because the report was rendered unconditionally — making the failure path
  honest is what exposed it.
- **The selection column was blank.** The marks were `[x]`, `[+]` and `[ ]`,
  and a `DataTable` cell given a `str` is parsed as Rich markup: `[x]` is a tag,
  not text. An interface for choosing tools showed no choices.

The lesson is the one this file keeps relearning, in a new place: **a second
implementation of a path that already works is where the bugs go.** The unit
tests were no help — 92 of them passed throughout, because not one of them
called `execute()`. There are now 100, and all 8 of the new ones fail against
the pre-fix source (checked by stashing `src/` and re-running them).

Two of the eight need real Snakemake and are skipped where it is absent, which
is deliberate: the defect was in Snakemake's own field names, and a fake event
stream would have agreed with whatever the code assumed.

### The default database directory was `./databases`, so it was per-run
A cwd-relative default meant the same command run from two directories fetched
a second copy of a set whose measured part is 143.2 GB. Nothing warned, and
nothing pointed the two runs at each other: the flag existed, so the fix on
thylakoid was to type `-d /evo/postdoc/cm2-databases` on every invocation and
write in STATUS.md that forgetting it re-downloads 6.9 GB. A default that has to
be overridden by hand on every run to avoid wasting 143 GB is the wrong default.

Now `~/.comparem2/databases`, with precedence `-d`, then
`$COMPAREM2_DATABASES`, then that. Home-relative rather than under `--output`
because databases outlive any one run's results, and outside any checkout so
deleting a checkout does not cost a re-download — which is the rationale
STATUS.md had already recorded for placing them by hand.

`$COMPAREM2_DATABASES` is a deliberate third way to say the same thing, against
the usual preference for fewer. It is the only one that can be set *once*: `-d`
has to be retyped, and a home directory is the wrong place for 143 GB on a
cluster with a quota on it. Without it the shared default is unusable exactly
where the sharing matters most, and the answer would be back to typing `-d`.

The resolved path is now printed with the download list, because the default is
no longer somewhere the user can see:

```
to download: checkm2, gtdb, bakta-light, amrfinder (143.2 GB + 2 of unknown size) -> /evo/postdoc/cm2-databases
```

Found alongside it: `render_report` was passed the *unresolved* `args.databases`
while `prepare` and the TUI got the resolved one. Harmless while the default was
relative to a cwd that had not changed, and a real divergence as soon as it was
not. A test now asserts the two agree.


### Nine minutes of CarveMe was a presolver eating the optimum
Carl asked whether CarveMe could be rewritten to run faster. It could not: the
Python is 8 s of 609. What the profile did show is that the 601 s in the MILP is
CarveMe's own hardcoded 600 s limit expiring, so the run was not solving, it was
giving up — and returning the best point it had found.

Then the same problem solved in 9.8 s, to proven optimality, under a different
SCIP build on the same machine. Byte-identical problem: the MILP was written out
from each build and the md5s matched (`fba2ad10…`), same carveme 1.6.6, same
DIAMOND hits.

| SCIP | MILP | status | model | annotated rx dropped |
| ---- | ---: | ------ | ----- | -------------------: |
| conda-forge 10.0.3, as shipped | 601 s | timelimit | 1,193 rx / 750 met | 253 |
| the same, `presolving/milp/maxrounds=0` | 42.7 s | gaplimit | 1,742 rx / 1,175 met | 44 |
| PyPI wheel 10.0.2 | 9.8 s | optimal | 1,738 rx / 1,176 met | 45 |

The one structural difference between the builds is PaPILO, which conda-forge
links and the wheel does not. Two hypotheses tested and rejected on the way: not
the SCIP version (conda-forge 10.0.2 has PaPILO 3.0.0 and was still in the
MILP when it was stopped at 300 s — it was never run to a conclusion, so "just
as slow" is an inference from that, not a measured end state), not symmetry
handling (`misc/usesymmetry=0` ran to the 900 s limit).

**The shipped run is also worse, not just later** — 253 of 1,069 annotated
reactions dropped against 45 — so "the models were fine, just slow" was never
available.

Reduced to the `scip` command line, one build (conda-forge 10.0.3 / SoPlex
8.0.3 / PaPILO 3.0.1), one written-out problem, `limits/gap 0.001`:

| run | reported | time | primal | dual |
| --- | -------- | ---: | -----: | ---: |
| `read lp; optimize` | time limit | 300.0 s | 913.500 | 935.696 |
| `read lp; read sol_947; optimize` | gap limit | 4.3 s | 947.500 accepted as feasible | 947.796 |
| `read lp; read sol_943; optimize` | **optimal** | 5.1 s | 943.500 | 943.500 |
| `read lp; set presolving milp maxrounds 0; optimize` | **optimal** | 7.5 s | 947.500 | 947.500 |

**This corrects the first reading of the same measurements**, which took the
943.4997 the PyPI wheel called optimal to be the optimum and PaPILO to be
straightforwardly unsound. It is not that clean. 947.4997 is feasible too, so
the wheel's "optimal" was also wrong — and **rounding either point's binaries to
exact integers leaves an infeasible LP**, checked independently of the solver's
own checker. Both are feasible only at tolerance: 5e-11 on the constraints,
exactly 1e-6 on integrality, against a default `feastol` of 1e-6.

So the honest statement is that "optimal" is not well defined on this problem at
default tolerances, because of CarveMe's conditioning: `minmax_reduction`
couples every flux to its indicator with bigM=1e3 against eps=1e-3, six orders
of magnitude apart. A SCIP maintainer would be within their rights to call the
dual bound a scaling consequence rather than a bug, and the report upstream is
scoped accordingly — the reproducible, tolerance-independent claim is that one
build and one file yield three different objective values, two of them labelled
optimal, and that this presolver costs 12–60x wall time and a quarter of the
annotated reactions.

Both genomes tested behave the same way — E8202, 3,185 proteins: 908.1 s and 261
annotated reactions dropped, against 16.7 s and 45. Both are *E. faecium*, which
is the limit of what has been checked.

Worth reporting upstream: any conda install of CarveMe has this, and the models
it produces are quietly smaller than the method asks for. The 4 MB models this
repository recorded as CarveMe's output on 2026-09-02 are those models.

### `carve` was silently overwriting Bakta's feature table
Found while tracing the above. `carve` derives its DIAMOND output path from its
*input* path — `os.path.splitext(inputfile)[0] + '.tsv'`, in `carve.py` — so
`carve …/bakta/<sample>.faa` writes `…/bakta/<sample>.tsv`. Bakta writes a
feature table under exactly that name.

In the run of 2026-09-02 morning, `samples/116_2/bakta/116_2.tsv` holds
12-column DIAMOND hits against BiGG gene ids. Nothing caught it because nothing
depends on that file: Bakta declares only its GFF3 and its FAA, and the report
reads the GFF3. So a declared-outputs discipline stopped this from breaking a
run, and did nothing to stop it happening — the file was wrong for anyone who
opened it.

The wrapper links the FAA into CarveMe's directory and points `carve` at the
link, which puts the hits in `carveme/<sample>.tsv` and makes them a declared
output rather than a stray. `Tool.files` could not have done it: it maps a path
to *content*, and this content is another rule's output, which does not exist
when the Snakefile is rendered.

### The content-addressed environment sharing was executed, not just argued
`catalogue.py` had a paragraph explaining why AMRFinder's fetch rule and its
analysis rules end up in the same deployed environment — necessary, because
`amrfinder -u -d <dir>` refuses any directory but the default, so the database
can only live inside `$CONDA_PREFIX`. The argument was read out of Snakemake's
`conda.py` and never run: the 09-02 `--use-conda` run stopped at
`--until seqkit mashtree treecluster skani checkm2`.

Run today. **Measured**: 6 → 8 environments in the shared prefix, five rules
(one fetch, four analysis) into one directory
`68e8563502afe8a0983c6c2bb5b459c1_`, because the two rendered env files are
byte-identical at md5 `cb5de824e5b0359eeb580a51570bd742`. The fetch rule built
the environment and the analysis rules joined it, so the sharing does not
depend on ordering — which was the specific way it might have "worked by
accident". Output byte-identical to the pixi run for all four genomes.

The evidence that matters is not exit 0. The pixi environment also carries
amrfinder 4.2.7 *and* database 2026-08-07.1, so correct output would not
distinguish the deployed binary from pixi's. AMRFinder prints its own software
and database directories on every run, and the per-sample log names the
deployed environment for both. Tools that report their own paths are worth more
than tools that do not, when the question is which copy ran.

### The stale thylakoid checkout was stashed, not discarded
STATUS.md had recorded the exact commands to bring `/evo/postdoc/CompareM2`
current and had verified nothing in its rsync snapshot was unique. The
amrfinder run forced the issue — it needs master's `catalogue.py`, and the main
checkout is the only one with a `.pixi`, so the other clones were not an option.

Carl's call was `git stash push -- src tests` plus a `.rsync-snapshot-backup/`
directory for the three untracked files, rather than the `git checkout --` and
`rm -f` the earlier plan proposed. Same end state, recoverable. The general
form: when a cleanup on a remote machine is believed safe but the belief has
not been tested by needing the files back, the reversible route costs nothing.

### `--use-conda` needs conda on `PATH`, and on thylakoid that is not obvious
The first attempt died at DAG construction with `Error running conda info. Is
conda installed and accessible?`. `conda` is not in `.bashrc`, not in the pixi
environment, and not in `pixi.toml` — it is a pixi **global** tool at
`~/.pixi/bin/conda`. The 09-02 script's `export PATH=$HOME/.pixi/bin:$PATH` was
load-bearing for a reason that had nothing to do with finding `pixi` itself.

Not a code defect: a bioconda install has conda by definition. Recorded because
it costs a failed run to rediscover, and the error message points at the
machine rather than at `PATH`.

### Snakemake installs the tools, and the flag that said so is gone
`--use-conda` was a flag, defaulting off, and the pixi environment carried all
thirteen co-solvable tools. Two models: pixi for development and HPC, per-rule
conda deployment for distribution.

Carl's call, and it is the right one: **pixi is how you install a development
environment, not a deployment model.** Snakemake can install the tools, so it
should, always. Deleted rather than automated — an earlier proposal in the same
conversation was a tri-state `--use-conda auto` that would guess from whether
the tools were on PATH, which is a heuristic standing in for a choice that
should not exist.

What went with the flag: `--isolated-launcher`, `Tool.isolated`,
`Tool.executable`, the `per_rule_conda` and `launcher` parameters threaded
through three modules, `missing_executables()`, and the thirteen tool
dependencies in `pixi.toml`. 417 insertions against 403 deletions across eleven
files, and the net effect is that a user types `cm2 *.fna`.

Two things this fixed that were not the point:

- **The tool set had been pinned in two files, and they had drifted.**
  `gtdbtk>=2.7` in `catalogue.py`, `gtdbtk = "*"` in `pixi.toml` — and the
  `osx-arm64` probe recorded in STATUS.md had already shown what the unpinned
  spec resolves to: gtdbtk 1.0.2 from 2019. `catalogue.py` is now the only
  place a tool is named.
- **`CheckM2` stopped being special.** `isolated=True` existed to give one tool
  a `conda:` directive when the others had none. Every rule has one now, so the
  DIAMOND 2.1.x-against-2.2.x conflict needs no mechanism — only a second
  environment, and a comment saying why.

The preflight shrank from fourteen checks to one. Nothing is expected on PATH,
so what is left to check is `conda` itself — which is the failure that actually
fires, and did, an hour earlier the same day.

### Two environments, not fourteen
Content addressing makes an environment per tool nearly free to write, and that
is the trap: fourteen solves cost fourteen copies of DIAMOND, python and numpy,
paid on first run. Eight single-tool environments had already measured 8.6 GB.
It is v2's 25 environments in a cheaper disguise.

So environments are **named** rather than derived per rule. `Tool.environment`
and `Database.environment` select a file; `Tool.conda` is the whole package list
of that environment, not the tool's own package. Eighteen rules, two names.
`render_envs` raises if one name is ever given two different package lists,
because that would write one file twice and let whichever rule rendered last
decide what the other one ran in.

**Measured 2026-09-03**, fresh prefix, all fourteen tools: 7.7 GB — `main`
6.0 GB, `checkm2` 1.8 GB — against 8.6 GB for eight of them separately. Both
solves took 76 s with a warm package cache, and 29 of the 30 rule activations
went to `main`.

Two rather than one is the DIAMOND conflict and nothing else. Verified
2026-09-01: `bakta>=1.10` co-solves with all thirteen others and fails only when
checkm2 joins.

**The thirteen-way co-solve was not a new risk.** It is the same solve
`pixi.toml` carried and every verification run used — 422 packages, seqkit
2.13.0, bakta 1.12.1, panaroo 1.8.0, gtdbtk 2.7.2, DIAMOND 2.2.5. What it does
change is the blast radius of an unconstrained spec: in a thirteen-way solve any
one tool can be reached backwards to satisfy another. So **every tool now
carries a `>=` floor**, at the build verified on linux-64, and a test enforces
it. Three had one before. `mlst` is floored at 2.33 rather than its verified
2.35.0, because 2.34+ need the Linux-only `libxcrypt1` and the lower floor keeps
the macOS finding in STATUS.md reproducible.

### The seven tools that had never been deployed
Making conda deployment the only model meant the seven tools never executed that
way had no verified route at all: gtdbtk, mlst, panaroo, snp-dists, fasttree,
carveme, biosynthesis. Two were the ones worth worrying about, because both
depend on *which interpreter a rule's shell gets* — `gtdbtk`'s post step runs
our own code through an absolute `sys.executable`, and `carve_scip.py` runs
under a bare `python` from the tool's environment, the mirror image. Both were
reasoned about in `catalogue.py` for exactly this case and neither had been run
in it.

Executed the same day: 31 of 31 steps, exit 0, all fourteen sections in the
report. Results agree with the pixi run — byte-identical for seqkit, amrfinder,
checkm2, gtdbtk, mashtree and treecluster; identical but for an embedded
absolute path for mlst and skani; identical values in a different row order for
snp-dists, whose order comes from panaroo's non-deterministic alignment; same
topology with sixth-decimal branch lengths for fasttree; same 3,780 clusters and
2,091 core for panaroo. Numbers in STATUS.md.

### `cm2 --setup`, and the install-time version that is not possible
The first run pays for the environment build, and it pays *before the first
job* — Snakemake deploys during DAG construction, so the run is silent for the
duration and looks hung. Carl asked what it would take for the installation to
have done it already.

**It cannot be the installation, and the reason is measured rather than
policy.** Snakemake addresses a deployed environment by
`md5(realpath(conda_prefix) + env file content)`, so the prefix is part of the
environment's identity. One byte-identical `main.yaml` went to three prefixes
and got three directories: `f35bbb1f…`, `00dbdb48…`, `efd5ffa1…`. The prefix is
a runtime choice — `--conda-prefix`, or `$COMPAREM2_CONDA_PREFIX`, which on a
cluster is set by everyone because home has a quota. A package deploying into
the default location at install time would therefore have built the wrong
directory for precisely the users who need it most.

Two further reasons, **not verified here** and flagged as such: `post-link.sh`
is the only hook and it runs inside the conda transaction that is installing
CompareM2, so it would be driving conda recursively against the same package
cache; and bioconda discourages post-link scripts, expecting them fast,
offline-safe and prefix-local. A 7.5 GB, 62 s post-link is none of those.

So it is an explicit step instead, and Snakemake already had the mode:
`--conda-create-envs-only`. Two properties make it usable as *setup* rather
than as a first run in disguise, both measured before the code was written:

- **No assemblies.** The DAG needs the per-sample FASTA at the bottom of it to
  exist, so a 26-byte stub goes into a temp directory and is deleted. Nothing
  reads it, because nothing runs.
- **No databases.** `-d` pointed at a directory that does not exist and the DAG
  still closed, because every database path in it is some rule's output. Setup
  works before any of the 62.5 GB has been fetched — which is the whole point,
  since otherwise "set up first" would mean "download 60.8 GB of GTDB first".

The temp workdir is safe because the *output* directory is not in the hash —
also measured: two runs with different `--output` against one prefix, and the
second built nothing.

Numbers: 61.7 s cold for both environments, 1.97 s when they exist, zero tool
outputs, no scratch left behind. The acceptance test is that a later real run
reuses it, and a `--until seqkit skani` run activated the same
`efd5ffa1fbfe3b0c3288c33676aaf20a_` and finished in 3.5 s.

`inputs` had to become `nargs="*"`, which means argparse no longer catches a
bare `cm2` — so that case now raises its own message naming `--setup`.

### Databases to /midifiler, environments to /evo
Carl's call, correcting a proposal of mine that put both on `/evo`: "the evo
drive is not very large." Measured — `/evo` is NVMe with 785 G free, and the
database root is 101 GB on disk. `/midifiler` is 13 T with 3.0 T free.

So the two variables now point at different volumes on purpose:

    COMPAREM2_DATABASES=/midifiler/carl/cm2_db_v3      13 T, spinning
    COMPAREM2_CONDA_PREFIX=/evo/postdoc/cm2-envs-two   NVMe

Databases are 94 GB of GTDB read sequentially, which is what a spinning volume
is good at. The environments stay on NVMe for three reasons: 7.5 GB does not
help the space problem, a conda environment is tens of thousands of small files
and every rule activates one, and **moving the prefix invalidates every
environment** because `realpath(conda_prefix)` is in Snakemake's hash. The
database path is in no hash at all, which is what made the move safe — verified
the same day, when `--setup` against a nonexistent `-d` and a real run against
the true one shared an environment.

Moved with `rsync -a` rather than `mv`: cross-filesystem, so `mv` is a
copy-then-unlink with no resumability and no verification. 101 GB in 12m23s at
138 MB/s, then `rsync -ani --delete` returned nothing and both roots measured
107,812,346,055 bytes. `bakta` is a *relative* symlink to `bakta_dl/db-light`
and survived. Proven working before anything was deleted: a full dry-run listed
only `amrfinder` to fetch — its marker is per-run by design — and a real
checkm2 run completed against the 2.9 GB database at the new path.

### `.bashrc` was pointing v3 at v2's databases
Found while checking whether the defaults were reasonable. Three stacked
`COMPAREM2_DATABASES` exports, v2-era, last one winning:
`/midifiler/carl/comparem2_databases`. It exists, holds 480 GB of v2 data
(`cm2_v2.15`, `cm2_v2.16`, `gtdb_sketch_release226`), and carries **none** of
v3's ready markers — no `bakta/version.json`, no `checkm2/checkm2.dmnd`, no
`gtdb/metadata/metadata.txt`.

So an interactive `cm2 <genomes>` with no `-d` would have seen all four
databases as absent and started refetching **62.5 GB**, including GTDB's
60.8 GB, into v2's directory. `/midifiler` has 3.0 T free, so it would have
succeeded — silently, slowly, and in the wrong place.

It never bit because every v3 verification run passed `-d` explicitly. That is
the shape of the problem: a default that is only ever exercised by people who
are not testing. The v2 lines are commented rather than deleted, since the
variable name is shared and a v2 run may still want them; `COMPAREM2_PROFILE`
is left alone because v3 never reads it. Backup at `~/.bashrc.bak-2026-09-03`.

**The defaults themselves are the right shape and the wrong place.** Both fall
back to `~/.comparem2/`, which is shared across runs and outside any checkout —
that shape was a fix, since the database default was once `./databases` and two
runs from different directories fetched up to 143 GB twice. But home here has
101 G free against ~108 GB of databases and environments, so the defaults would
just fit and be a bad idea; on a cluster with a quota they fail partway through
a 60.8 GB download. Which is the whole reason both variables exist.

---

## 2026-09-04

### `carve_scip.py` prints SCIP's status, gap and bounds for every solve
Instrumentation, not a fix. The 2026-09-02 seven-genome *S. aureus* run
(measured) had three genomes come back with 1,236–1,263 reactions against
1,609–1,636 for the other four — 22–24% short — at 613–619 s against 22–30 s,
and nothing said so. Reading the source afterwards explains why nothing could:

- `carveme/reconstruction/carving.py:180-184` sets `limits/time=600` **and**
  `limits/gap=0.001`, then calls `solve(allow_suboptimal=True)`. The pair is
  the design: stop at a 0.1% gap, or at 600 s if the gap has not closed. It
  applies to SCIP only — the Gurobi and CPLEX paths get no ceiling.
- `reframed/solvers/scip_solver.py:15,20` maps `timelimit` and `gaplimit` alike
  to `Status.SUBOPTIMAL`, and `solver.py:137` accepts SUBOPTIMAL when
  `allow_suboptimal` is set. So a model stopped by the clock is returned through
  the same path as one that converged, with no exception and exit 0, and the
  status CarveMe can see does not distinguish them.

The existing `SCIPSolver.solve` patch is the only place that still holds the
SCIP problem after a solve, so the report goes there. It prints and does not
act: **which of the two stopping criteria to move is not yet answerable.** If
the gap at cut-off is small then `limits/gap=0.001` is the wrong knob to be
strict about on this problem — the wrapper's docstring already argues that
"optimal" here describes solver tolerances rather than the network — and
relaxing it costs no packaging change. If the gap is large, the incumbent is
genuinely poor and the PyPI SCIP wheel (measured: 10 s and `optimal` on the
*E. faecium* instance where conda-forge gave 601 s and `timelimit`) is the
candidate, at the cost of a `pip:` entry in a conda environment. Nobody has the
gap number yet, which is the point of the line.

Raising the ceiling is already ruled out: N315 re-solved at 3,600 s still
returned `timelimit`, with 7 reactions *fewer* than the 600 s run (measured).
A commercial solver is not a candidate for a pipeline other people install.

Still undone, and separate: a truncated model reaches the report and
`biosynthesis` indistinguishable from a converged one. The status is not yet a
declared output, so nothing downstream can refuse or flag it.

### The CarveMe timeout is the instance, not the solver build — both candidates are dead
Measured 2026-09-04, numbers in STATUS.md. The entry above proposed two
branches, chosen by the gap at cut-off. The gap turned out to be 1.0–2.1%, and
**both branches are ruled out**:

- Relaxing `limits/gap` fails for the reason the small gap suggested it would
  work. At 2% the solver would stop *earlier*, on an incumbent just as sparse.
  The gap being small is not evidence that the answer is nearly right.
- The PyPI SCIP wheel fails outright: the same three genomes hit the same 600 s
  ceiling with the same ~1,245 reactions, PaPILO absent. The presolver fix was
  real on 116_2 and does not generalise, so **the slow instances are a property
  of the genome, not of the build.**

The number that reframes it: MRSA252's 1,251-reaction incumbent scores 920.4
where NCTC8325's *converged* model scores 907.1 with 1,629 reactions. CarveMe's
objective does not reward reaction count, so a third of the network can come or
go inside 1.5% of objective. Which model you get depends on the search path.
That is the formulation — bigM 1e3 against eps 1e-3, integrality feasible only
at `feastol` — and not something a solver swap or a longer clock addresses.

**So this stops being a fix to find and becomes a limitation to disclose**, and
the disclosure work that was filed as secondary is now the whole of it: the
solve status has to become a declared output, the report has to show it, and
`biosynthesis` must not read a truncated model as though it were converged.
Reversal is deliberate — the earlier entry named the wheel as the likely
answer and it is measured not to be.

One measurement would still change the severity, and has never been made:
whether the sparse models drop reactions the annotation supports. That column
is what convicted the shipped build on 116_2 (253 against 44) and is unknown
for these three. Small is not the same as wrong.

### The release after 3.0.0 is 3.1.0, not 3.0.1
Carl's call, 2026-09-04, having asked for 3.0.1 and been shown the objection.
`cm2 --demo` is a new user-facing flag, and the README, quick start and usage
pages are rewritten around it — a patch number would say "nothing new here"
about the one thing the release exists to ship. Nothing downstream is affected
either way: the recipe's `run_exports` pins at `max_pin="x"`, so 3.0.1 and
3.1.0 pin identically.

The other half of the release changed with it. **Pushing the tag is now the
whole release**: bioconda's autobump is on, and unlike the v2→v3 bump it is the
right tool for this one, because the published recipe is already the v3 shape
and a version-and-checksum bump is the entire change. The bot's #68821 was a
hazard only because it was generated from v2's recipe while the shape change
was in flight; it is closed. So no hand-written PR this time — see
[recipe/README.md](recipe/README.md), where the release steps now split at that
line.

What was tested first is in [STATUS.md](STATUS.md). A macOS check cannot run a
tool and never will, because `Tool.conda` renders the whole thirteen-tool
environment for any subset — so the laptop covers everything up to the first
job (DAG, extraction, wheel, entry points, both docs checks) and nothing after
it.

**Amended the same evening:** that was first written as an argument — the
tagged commit is two prose strings from the tree that ran 11 of 11, so its
tools "need not" be re-run. Carl asked why it had not simply been run, which
was the right question, and the answer is that ssh to thylakoid cost four
seconds: 11 of 11 at `v3.1.0` detached, identical seqkit md5 and identical
skani, mashtree and treecluster numbers. An argument that a run is unnecessary
is worth less than the run whenever the run is cheap, and here it was cheaper
than the paragraph defending it.

---

## 2026-09-07

### `--profile`, and the TUI submits too
The pipeline could not submit a job. Both execution paths were pinned to local
execution — `cli.py` built a Snakemake command with no `--executor`, and
`runner.py` passed `executor="local"` literally — while the docs claimed
`--cores 64` would submit. It starts 64 processes on the machine you typed it
on.

`--profile` is a passthrough. Cluster submission belongs to Snakemake and
already worked there; both executor plugins had been declared dependencies all
along, so SLURM, PBS, SGE and LSF arrive with the flag and no new code. v2
spelled the same thing `$COMPAREM2_PROFILE` and shipped fourteen profiles
in-tree — two of the four it advertised as cluster-specific were byte-identical
to the templates, placeholder account and all — so v3 ships none and documents
one.

The first proposal was to *refuse* `--profile` with `--tui`, on the grounds
that the API has no notion of a profile. Carl rejected it: watching a queue
from a frontend is exactly when a progress display earns its keep. The profile
branch therefore calls `snakemake.cli.parse_args()`/`args_to_api()` — Snakemake's
own CLI, in-process — because `config.yaml` is read as argparse *defaults*, and
hand-mapping that onto `execute_workflow(executor=..., executor_settings=...)`
would be re-implementing a parser and getting a subset right.

**A measurement retired the reason for one of the guards.** Withholding the
default `--cores 4` under `--profile` was justified in the code as preventing a
cluster run being capped at four jobs. Measured on GenomeDK: under `--cores 1`
against a profile's `jobs: 20`, all three test jobs still started at 2.2 s.
`--cores` governs local scheduling, not submission. The behaviour is kept for
the smaller true reason — a profile may set `cores:` itself and our default
would replace it — and the comment now says the measured thing.

### Downloads are `localrules`, because Snakemake decides where a rule runs
Carl's observation that running `--setup` on a frontend is obvious was correct,
and chasing it found the thing that is not: the four `download_*` rules are
*rules*, so under a profile Snakemake submitted them like everything else.
GenomeDK's compute nodes have no outbound network, so the 60.8 GB GTDB fetch
would have been sent to the one machine that cannot reach the internet. No
amount of user discipline fixes that.

Verified: given `localrules: fetch`, `fetch` ran on `fe-open-01` while a sibling
rule went to `cn-1050` in the same run.

### Two environments became six, because two stopped solving
**This reverses "Two conda environments, and adding a third needs a reason"**,
which stood from 2026-09-01. The reason turned out to be the strongest kind:
the thirteen-tool `main` environment stopped solving at all. Not on one machine
and not from new configuration — the identical spec built a working 6.0 GB
environment on thylakoid on 09-03 and, re-solved there on 09-07, failed after
6 min 46 s; on GenomeDK after 4 min 07 s. `perl-bioperl` could not be placed,
taking `mlst` and `panaroo` (through `prokka`) with it, while `mashtree` wanted
`perl >=5.32.1` and `perl-bio-samtools` was offered only as perl 5.26/5.22
builds.

Three things were ruled out by measurement before the split, and are recorded
so they are not re-tried: the `_python_rc` / python-3.14rc line is a red
herring — an earlier note blamed it and that was a misreading of the solver's
own tree; pinning `curl` and `tar`, the two specs with no floor, fails
identically; `channel_priority: strict` fails identically.

The lesson is about co-solving, not about perl. Thirteen tools in one
environment means thirteen sets of transitive constraints that must hold at the
same moment, so one ecosystem going bad upstream takes the other twelve down.
Every candidate group solves in under 30 s where the thirteen-way solve fails
after four minutes, so the split is by *ecosystem*: `basic`, `perl`,
`annotation`, `gtdbtk`, `carveme`, `checkm2`.

An environment per tool is still wrong — v2's 25 in another form — and the
floors are still mandatory for the same reason as before.

Two things fell out of it that were not the goal. A subset now builds only what
it needs: `--until seqkit skani` builds one 55 MB environment where it used to
build the whole thirteen-tool `main`. And all six together are 1.4 GB by `du`,
against 7.7 GB recorded for the old two — but conda hardlinks shared packages
and the two figures were measured differently on different machines, so that is
not a like-for-like comparison and is not claimed as one.

`envs/locks/` holds `conda list --explicit` from thylakoid's surviving
environments, 394 packages and 130. Nothing reads them. They were captured when
the split was still one of two options, and they stay as the record of a set
that ran — the lock installs on GenomeDK in 6.9 s where the floors-only solve
failed in 4 min 07 s, which is worth knowing if drift ever hits a group that
cannot be split further.

### The TUI reads the output directory before it shows anything
Opening `--tui` on a directory that had already been run in showed fourteen
rows of `pending`. The state dict was seeded `{t.name: PENDING for t in
CATALOGUE}` and nothing ever looked at the filesystem, so the interface whose
job is to say what still needs doing could not answer that question at all —
the only way to find out was to press `r` and watch Snakemake skip things.

The signal is the tools' **declared outputs**, which is the same thing
Snakemake's resumability is decided on, so what the table says on startup is
what a re-run would actually skip. Three places were already asking that
question in three different spellings — `cli.any_outputs_exist`,
`report.render_report`, and nothing at all in the TUI — so it became one
function, `tools.completion()`, counted per unit of work because the three
callers need three different answers from the same files:

| caller | question | reads |
| --- | --- | --- |
| TUI | which analyses are done | `units == complete` |
| `any_outputs_exist` | is there anything to report | `complete > 0` — one finished genome counts |
| report | does this section render | `started > 0` — partial runs stay readable |

`already run` is a separate state from `done`, not an alias: "done" is
something the user watched happen this session, and putting this session's name
on a file left by a run last week is the kind of small false claim that gets
believed. `part-finished` is its own state too, and the startup line says those
will be redone — missing one declared output is exactly what makes Snakemake
re-run a rule.

Two consequences fell out that were not the goal. A tool that is *not* selected
still reports results it has on disk, because the mark column already carries
the selection and `not selected` over a finished analysis reads as "there is
nothing there". And pressing `r` in a finished directory used to print
**"Nothing ran. No report written."** over a complete set of outputs: Snakemake
emits no job events for a rule it skips, so every row settled to `not run` and
the report was withheld. The TUI now calls `any_outputs_exist` for that
decision — the CLI's own answer since 09-02 — rather than counting the events
it happened to see, so the two paths cannot disagree about whether there is
anything to report.

### The TUI shows where the databases and environments come from
`$COMPAREM2_DATABASES` and `$COMPAREM2_CONDA_PREFIX` are exported once in a
shell profile and then never looked at again, and neither appeared anywhere in
the interface. The failure mode is on record above: three stacked
`COMPAREM2_DATABASES` exports in `.bashrc` had v3 pointed at v2's database
directory on 2026-09-03, which would have silently refetched 62.5 GB into the
wrong place. Neither variable produces an error when it is wrong — one costs a
re-download, the other a re-solve of every tool environment, because Snakemake
keys a deployed environment on the prefix's realpath.

So four lines above the tool table: output, databases, tool envs, execution,
each with **where the value came from** rather than just the value.
`run_settings()` lives in `cli.py`, next to the defaults that read those
variables. The case worth naming is the third one — a variable that is set *and*
overridden by a flag looks identical to a variable that was never set, so that
reads `given, overriding $COMPAREM2_DATABASES`. `execution` says `local` or the
profile, since "am I actually submitting to the queue" is the same class of
question.

### `--unlock` takes no assemblies
`pixi run comparem2 --unlock` exited on `no assemblies given — pass one or more
FASTA files`, which is a demand for input that clearing a lock never reads.
The flag was handled at the bottom of `main`, after the input requirement,
after canonicalisation and after `prepare` — so the one working invocation,
`comparem2 *.fna --unlock`, copied every genome into the workdir and re-rendered
the Snakefile before releasing anything. On a killed 60.8 GB download that is
the wrong order of operations twice over.

A lock belongs to an output directory. `--unlock` now returns early on
`--output` alone, reusing the Snakefile the dead run already left in
`<output>/.comparem2/` — if there is none, no run ever started there, so there
is no lock, and it says that instead of handing back a Snakemake traceback.
Assemblies are still accepted and still ignored: adding `--unlock` to the
command that just died is how anyone reaches for it, and unlike `--setup` and
`--demo` there is no risk of the command looking as though it analysed them.

### The TUI animates while it runs
"When you start the run there is no way to see if the screen is frozen" —
Carl, and correct: between pressing `r` and the first `job_started` event
nothing on screen changed at all. That gap is the *longest* one in a run,
because a first run solves six conda environments before Snakemake emits a
single job event, and the interface deliberately quietens Snakemake's own
"Creating conda environment" output because it scribbles over the display.

One line above the progress bar, at 10 fps: a spinner frame, the tools that are
running, and the elapsed clock. Two decisions in it worth keeping:

- **Driven by a `set_interval` timer on the UI thread, not by the worker.** That
  is what makes it evidence rather than decoration — if the interface is
  genuinely blocked, the spinner freezes with it. A spinner animated from the
  Snakemake thread would keep turning through exactly the failure it is there
  to detect.
- **When nothing is running the line says which nothing it is.** `starting up —
  the DAG, and tool environments on a first run` before the first job, `no job
  running — waiting on Snakemake` after it, `collecting outputs and writing the
  report` once the event stream has ended. Rendering the report reads every
  output and takes seconds on a real run, so without the third the line would
  have been claiming a conda solve at the moment a user is most likely reading
  it.

The progress bar goes indeterminate — it pulses — from `r` until the first
`progress` event carries a total, and reverts to a drawable total when the run
ends. An indeterminate bar animates for as long as it has no total, so leaving
one going after the run would be the same lie in the other direction. A run
that reported real numbers keeps them: a failure at 1 of 4 stays at 1 of 4
rather than being reset.

At the end the line is replaced by `not running — the last run took 4m 12s`. A
stopped spinner and a hung interface look identical, so the animation is
removed rather than frozen mid-frame.

Six tests, and the two that matter are about the honesty rather than the
motion: that the animation is torn down on every exit path including one that
raises, and that a spinner frame is never also a status glyph — `◐` already
means `part-finished` in the column two lines above.

### The TUI opens with nothing selected
"I don't think that all tools should be enabled as default. Let the user decide
what to run manually" — Carl. It had seeded the selection with all fourteen
tools whenever `--until` was absent, which made `r` a 60.8 GB download for
anyone who pressed it before reading the second line, and made the interface a
confirmation step for a decision already taken rather than the place the
decision gets made. `a` still selects everything in one keypress, and `--until`
still seeds, so `--tui --until mashtree treecluster` opens on exactly those.

The change is one line; what it cost is the second one. `r` on an empty
selection had been a silent no-op, which was harmless when empty was a state
you had to press `n` to reach and is the worst available answer now that it is
the state the interface opens in — the key that runs things appearing to do
nothing at all. It now says what is missing and which key fixes it, and the log
says so on startup too.

The CLI default is untouched: no `--until` still means all fourteen. A command
naming no tools is unambiguous about wanting the lot, and there is no table in
front of the user to choose from.

### The TUI can unlock the directory it opens on
`--unlock` existed but only as a second command: quit the interface, remember
which `--output` the dead run used, retype it. The lock is also the one
condition that makes everything else the table says unreachable, and Snakemake
only reports it several seconds into a run that has already announced it is
starting up.

So the interface reads `<output>/.snakemake/locks/` on opening, says so if
anything is there, refuses `r` while it is, and `u` clears it. `cli.unlock()`
is shared with the flag — it returns the problem as text rather than raising
`SystemExit`, which is what the TUI can use, and captures Snakemake's output
because anything written to stderr scribbles over a Textual display.

**It asks first, and the dialog is about the one thing that cannot be known.**
A lock file holds a list of paths and no process id, so neither the user nor
this code can distinguish a killed run's lock from a live run's — and clearing
a live one puts two Snakemake processes on the same outputs. The dialog says
that in those words rather than asking "are you sure". Two further honesty
rules in the same feature: `u` during a run refuses, because that lock is this
run's, and after `snakemake --unlock` returns 0 the directory is re-read rather
than declared clear, since the exit code is not the same statement.

Detection is presence of `*.lock`, not Snakemake's own `Persistence.locked`,
which asks whether the locked paths intersect *this* DAG's files and therefore
needs a built DAG — unavailable before the run starts, which is exactly when
the question is being asked. Every run in one output directory is the same
workflow over the same outputs, so presence is the right answer here.

Verified against a real `snakemake --unlock` and a real generated Snakefile, on
a lock written by hand in Snakemake's format: `r` refused, `u` opened the
dialog, `n` left the lock, `y` removed it, and `--unlock` did the same from the
command line. **The lock was not one a killed Snakemake left behind** — see
[STATUS.md](STATUS.md).

## 2026-09-07 — the interface is legible on eight colours, and quitting says what it costs

### The selection was invisible over SSH, and a colour could not fix it
Reported from a real session on GenomeDK: with the command palette open, no row
looked selected. Not a rendering accident — arithmetic.

tmux ships `default-terminal screen`, an eight-colour TERM, so Rich renders
through its `standard` colour system and every RGB colour in the theme is
downgraded to one of eight. Textual's *blurred* cursor is `$primary` at 30%
alpha: `#0178D44C` blended over the surface is `#153854`, which downgrades to
ANSI 8, against a surface that downgrades to ANSI 0 — two near-blacks — and
`block-cursor-blurred-text-style` is `none`, so nothing else distinguishes it.
The palette is where it showed first because its list is `can_focus=False` and
is therefore *always* drawn blurred. Every built-in theme has this, including
the two ANSI ones.

**The fix is an attribute, not a colour: `text-style: reverse`.** A colour
cannot be chosen safely here — `ansi-dark` and `ansi-light` set `surface` to
`ansi_default`, which is whatever the user's terminal background happens to be,
so nothing can be guaranteed distinct from it. `reverse` is SGR 7: it inverts
whatever the row already is, in every colour system and every theme. Verified
at the byte level — a `standard`-system console emits `\x1b[7;36;40m`, reverse
plus cyan on black — and by the resolved styles in the real app, where the
palette's highlight comes out opaque `#1E1E1E`/`#0178D4`, ANSI 0 against ANSI
6, with `reverse` set. Focus is carried by `bold` on top, for the same reason a
second colour would not do.

Measured candidates before settling on this: solid `$primary` survives the
downgrade in `textual-dark` (6 vs 0), `textual-light` (8 vs 15), `nord` (7 vs 8)
and `gruvbox` (7 vs 8), but darkened variants collide with the surface in nord
and gruvbox, and no colour at all works against `ansi_default`.

### `q` asks, and says what happens to the jobs
Quitting mid-run was silent about the only thing that mattered. What it
actually did was measured with a probe of this module's shape (Textual 8.2.8):
`App.run()` returns in 1.52 s, the worker's next `call_from_thread` raises
`App is not running` about a second later — so Snakemake stops being driven
almost at once — and **the child process survived its parent** and finished its
work twenty seconds on. So: nothing downstream starts, no report is written,
the lock stays, and the jobs keep going.

The dialog now says which of two situations the user is in, because the answers
differ: with a profile the jobs are in a queue that does not care that the
frontend is gone, and without one they are child processes of it. `s` quits and
stops them.

**Cancellation cannot go through Snakemake, and the reason is structural.** Its
scheduler reaches `executor.cancel()` from one place — a `KeyboardInterrupt`
inside its own loop — and installs the SIGTERM handler that would get it there
inside a `try/except ValueError` that silently skips when the scheduler is not
on the main thread (`job_scheduler.py:170-175`). `runner.run()` puts Snakemake
on a worker thread, so that handler is never installed and there is nothing to
signal. Reaching for the executor object is no better: the profile branch goes
through `args_to_api()`, which returns a bool.

So `cancel.py` does each half directly, and they are not symmetric:

- **A queue** is cancelled with `scancel --name <run_uuid>`. The SLURM executor
  plugin submits every job of a run under the run UUID as its job name *in
  order to* make `--name`-based cancellation possible — its own comment says so
  — and announces it as `SLURM run ID: <uuid>`, which `_Capture` now reads. One
  call covers jobs submitted but not yet started, which is the case that
  matters most.
- **This machine** needs the process tree signalled by hand, because
  Snakemake's local executor `cancel()` is `self.pool.shutdown()`
  (`executors/local.py:260`) — it stops *scheduling* and waits for what is
  running. It does not kill anything.

Both run under a profile, not one or the other: the analyses are queue jobs,
but the four `download_*` rules are `localrules` and run wherever Snakemake is,
so a 60.8 GB GTDB fetch is a child process rather than a job.

Two smaller decisions inside that. The outcome is printed by `departure()`
*after* the app has closed rather than inside it — `scancel` is allowed a
minute by the plugin's own code and the SIGTERM grace is three seconds, neither
of which belongs on a screen the user has just asked to leave, and text in a
restored terminal can be scrolled back to. And the decision is recorded on the
app (`left_mid_run`, `stop_requested`) at the moment it is made rather than
read back afterwards, because `self.running` is cleared by the worker thread
that is still coming down as the app exits — reading it after `run()` returns
is a race.

`action_quit` is overridden rather than binding `q` to a new action, so that
everything which quits goes through the question: the key, the command
palette's Quit entry, and ctrl+c — which in Textual 8 does not quit but points
at whichever key runs the `quit` action.

**What is not verified:** `scancel` has never been run from this code. The id
capture is tested against a synthetic log record, and the local half against a
real process tree, but no queue has been cancelled — see
[STATUS.md](STATUS.md).

### The process walk reads `ps`, not `pgrep`, because a zombie is not a job to stop
The test above passed on macOS and **failed on all three CI Pythons** —
`1 failed, 247 passed` on Linux against 250 passed on the laptop, the first red
CI on this branch. The log alone identifies the cause: the pid it complains
about survived SIGTERM *and* SIGKILL, which only an already-dead process can.

A process that has exited holds its place in the process table until its
parent reaps it, and in that test the parent is asleep for 30 s and never does.
`pgrep -P` lists such a zombie on Linux and does not on macOS — that difference
is the whole of the platform split, and it was never about the signalling.

So the walk now takes one `ps -A -o pid=,ppid=,stat=` snapshot and drops
anything in state `Z`, rather than one `pgrep -P` per node. Three reasons in
that order: **`ps` reports the state and `pgrep` cannot**, so the check is only
possible this way without a psutil dependency; one snapshot cannot shift
underneath the walk where a call per node can; and it is one subprocess instead
of one per node. Zombies are still walked *through* in case a table lags,
though a dying process's children are reparented at once.

This was a real defect and not only a red test: `stop_local()` re-scans after
the grace period and SIGKILLs what is left, so on Linux **every** cancelled
local run would have reported the processes its SIGTERM had already stopped as
having "needed SIGKILL" — the one line the user sees after the interface
closes, wrong about the only thing it says.

The new test asserts the zombie is gone from `descendants()` *and* that `ps`
still lists it as `Z`, so it is the unreaped case being checked rather than a
table that moved on. Verified to bite on macOS too, by running it against a
copy of the module with the state check removed: `[7034]` where the fixed one
returns `[]`.

### The release after 3.1.0 is 3.2.0, not 3.1.1
Same question as 2026-09-04, same answer, and for the same reason: `--profile`
is a new user-facing flag and a patch number would say "nothing new here" about
the thing the release exists to ship. The rest of the diff argues it harder than
3.1.0's did — the deployment shape changed from two environments to six, the
TUI grew a run picker, an unlock, a quit dialog and `cancel.py`, and `cli.py`,
`snakefile.py` and `tools.py` all changed. Nothing downstream is affected either
way: the recipe's `run_exports` pins at `max_pin="x"`.

**Why it goes out the same evening rather than after 3.1.0 lands.** Bioconda
holds an autobump PR for three days from *its own creation* and then merges it
(the rule is in [recipe/README.md](../recipe/README.md), with the evidence).
v3.1.0's PR was opened 2026-09-04T20:15:57Z, so it becomes eligible at
20:15:57Z on 09-07, and the bot pushes a newer tag onto that same PR rather
than opening a new one. A tag pushed this evening therefore ships tonight;
waiting for 3.1.0 to merge first starts a fresh three days and lands 3.2.0
around 09-11.

That is a timing convenience. The reason it is worth having is that **3.1.0 is
a version that cannot deploy its tools** — it renders the two-environment
`main`, which stopped solving upstream on 09-07 — so the channel currently has
nothing that works past the first job, and 3.2.0 is the fix.

## 2026-09-07 — the download rules get rows, not a special case

### A database is a row, keyed by its rule name
The four `download_*` rules were rules Snakemake reported on and rows the TUI
did not have, and the failure was silent by construction: `mark()` looked up
`download_gtdb` in a `self.state` keyed by the fourteen tools, missed both its
lookups, and returned. Three readers of that dict were wrong at the same
moment — the table, the activity line (`no job running — waiting on Snakemake`
through a 60.8 GB fetch), and `action_quit`'s count, which produced **"no job
has started yet"** in the dialog whose entire job is to stop that sentence
being read as *nothing to lose*.

**Rows rather than a message.** Special-casing the wording would have fixed one
of the three and left the other two, and a fourth reader added later would have
been wrong again. The table already is the answer to "what is this run made
of"; a download was simply missing from it.

**Keyed by rule name, and the key is `Database.rule` in `tools.py`.** Not
`db.name`: `checkm2` names both a tool and a database, so the two would collide
on one key and the display text no longer identifies a row. Putting the string
on the spec rather than in `snakefile.py` is the load-bearing half — it is what
`runner.Event.rule` carries back, so it is not only the generator's business,
and two definitions drifting apart would be invisible from both sides. The
Snakefile would still run; the interface would just stop finding the row again.

That change made `_row_key()` wrong. It read the *Tool column's text* and used
it as identity, which held only while every row was a tool. It now reads the
row's own key.

### `space` on a database does nothing, deliberately
Nothing chooses to download GTDB — GTDB-Tk chooses it by needing it — so the
row carries the dependency mark `▨` and no checkbox. The guard is not cosmetic:
without it the name falls into `self.selected` and the next `closure()` raises
`KeyError: unknown tool: download_gtdb` inside a Textual worker. Confirmed by
removing the guard and watching the test fail that way.

### What a real terminal found that the tests had not
Driven under tmux against a synthetic event stream, which is the only kind
available on macOS. Two defects were visible in the pane and in neither the
unit suite nor the review that preceded it:

- an unneeded database read `to download`, a promise nothing in the DAG keeps,
  where the tool rows say `not selected` from the same position;
- `settle()` re-labelled all four databases at the end of a run, flipping an
  unneeded one back from `not needed` — undoing, one second later, the thing
  `sync_table` had been careful about.

Both now have tests. The general point is the one worth keeping: this interface
is where a wrong claim is *read*, and the unit suite checks the values behind
the claims rather than the sentence the user ends up looking at. A pane capture
is cheap and should be part of changing it.

### The cost line was fixed once already, on the other path
`cli.py` filters the download total by each database's `ready_path`, with a
comment recording the bug — "how `databases: 143.2 GB` came to be printed
before a run that downloaded nothing at all". The TUI's cost line never got
that filter and still totalled every database the selection needed. Carried
across, with the distinction the CLI does not have to make: `no databases`
means this selection needs none, and `none — all present` means it needs them
and they are here. Collapsing those two into one string would be the same
class of wrong claim as the rest of this entry.

`refresh_cost()` is now also called from `settle()`. It hung off `sync_table()`,
which only the selection keys reach, so the figure sat at 60.8 GB after the
very session that downloaded it.

---

## 2026-09-08

### A cluster profile is picked up from `$SNAKEMAKE_PROFILE`, and there is no `$COMPAREM2_PROFILE`
The question was how to make v3 pick up a profile from `~/.config/snakemake`
automatically, and whether that wants a `$COMPAREM2_PROFILE` like v2's. It
wants neither a new variable nor a search: **Snakemake already defines both the
default location and the variable**, and v3 was obeying one of them by accident.

`~/.config/snakemake/<name>/` is Snakemake's own default because
`appdirs.AppDirs("snakemake", "snakemake")` is what its profile lookup uses —
`get_profile_dir()` searches the cwd, that directory, and the system one, which
is why a bare `--profile slurm` already worked. And its `--profile` is declared
`env_var="SNAKEMAKE_PROFILE"` on a parser that subclasses
`configargparse.ArgumentParser`. Verified by running it, not only by reading it:
with `SNAKEMAKE_PROFILE` exported, `parse_args([])` comes back with
`executor='slurm'` and `cores=32` out of the profile's `config.yaml`
(snakemake 9.26.1; source read at the pinned 9.16.3 is identical).

**So the variable was half-honoured, in the worst available way.** The CLI path
shells out to Snakemake, which read the variable itself — so an exported
profile *did* submit the run, while `run_settings()` reported
`execution: local`, and a `--cores 4` nobody typed went out over the profile's
own `cores:` because that gate tested `args.profile` rather than the effective
profile. The TUI path uses Snakemake's API, which has no notion of a profile,
so the same variable submitted nothing at all. One exported variable, two
behaviours, and the header wrong in one of them.

The fix is to read it in `resolve_profile()`, the single place both paths get
their value from. A second name would have needed a precedence rule and bought
nothing; the honest version of "like v2" is to follow the upstream variable.

**`none` is the escape, and it is why every Snakemake we launch now names its
profile.** Snakemake spells "no profile at all" as `--profile none`, which
`resolve_profile()` turns into None. But *omitting* `--profile` on the
subprocess is not a local run — it is whatever the environment says. So
`profile_flag()` always emits one, `none` included, and `--setup` and
`--unlock` pin it to `none` outright: both are login-node bookkeeping, and a
profile would have changed the defaults they run under. Verified that this
costs nothing in the local case — with no variable set, `--profile none` and
passing nothing parse identically (`profile`, `executor`, `cores`, `jobs` all
None).

Nothing would have been submitted by those two either way:
`--conda-create-envs-only` reaches `dag_api.conda_create_envs()` and `--unlock`
reaches `dag_api.unlock()` in the same `elif` chain that ends at
`execute_workflow()`.

**What was rejected: scanning `~/.config/snakemake` and using what is there.**
A directory existing is not a decision to submit, and with `slurm` and
`slurm-gpu` side by side there is no right answer. Discovery is not worth a
surprise of that size.

**Unverified on a cluster.** 267 unit tests cover the resolution, the
attribution and every command line built; the GenomeDK submission of 2026-09-07
was through `--profile`, and no run has yet been started from the variable.

### A run Snakemake refused to start was reported as a run with nothing to do
Found by asking what `--profile path/to` does. The answer is three answers —
a directory with `config.yaml` submits, a directory without one and a path that
does not exist both hand off to Snakemake's own profile search, which fails
naming every directory it looked in — and the third of those is where the
defect was.

**Snakemake's CLI refuses a command line by calling `exit()`.** That raises
`SystemExit`, which is not an `Exception`, so `runner.run()`'s handler never
saw it and `threading` discarded it without a word. Measured:
`list(run(..., profile='empty'))` returned `[]` — no `done`, no `error`, no
job events, nothing.

The CLI path is unaffected: Snakemake is a subprocess there, the message goes
to stderr and the exit code is 1. The TUI path read that empty stream as
success. In a fresh output directory it recovered by accident and said
`Nothing ran`, without saying why. **In a directory holding earlier results it
said "every selected tool's output was already up to date" and wrote a
report** — about a run that had never started. `have` is satisfied by outputs
on disk and cannot distinguish the two, and neither can the table: "no job
events" means either *Snakemake skipped everything* or *Snakemake never got as
far as a job*.

Two changes, because the event alone was not enough. `runner.run()` catches
`SystemExit` separately and ahead of `Exception`, naming `--profile` as the
only part of that argv the user supplied — Snakemake's own message names the
directories it searched, but it goes to stderr, which is not where a TUI user
is looking. And the TUI keeps `run_error`, so the up-to-date sentence is
withheld and replaced by "The run did not start. The report below describes
what was already on disk."

Verified against real Snakemake, 9.26.1: the two failing profiles now yield
one `error` event each, and a valid profile still yields
`started, job_started, done`.

**What was not done: pre-validating the profile in `cli.py`.** Checking for a
`config.yaml` ourselves would duplicate a search that has already changed
between the pinned 9.16.3 ("no config.yaml found") and 9.26.1 ("no
profile.yaml (or config.yaml) found") and accepts `config.v<major>+.yaml`
patterns besides. `resolve_profile()`'s docstring already says to let
Snakemake's search produce that error; the bug was never the message, it was
that one execution path threw the message away.

### The variable was verified on GenomeDK the same day, and it submits
`SNAKEMAKE_PROFILE=slurm`, no flag, `comparem2 --demo`: 11 of 11 steps, 10 jobs
on `cn-1060` and `cn-1103`. The numbers are in [STATUS.md](STATUS.md).

Two things the local tests could not have shown. The bare-name form is what a
cluster actually uses — `~/.config/snakemake/slurm` was already sitting on
GenomeDK, written from the installation docs — and that path is the one
`resolve_profile()` deliberately does *not* resolve, so it had never been run
end to end. And the remote Snakemake is 9.16.3, the pinned version, where the
`env_var=` declaration had only been read rather than executed; the laptop's
venv is 9.26.1.

`--profile none` and `--setup` were checked in the same session with the
variable still exported: both stayed on the frontend and `sacct` recorded
nothing for either.

### The release is 3.3.0, and it rides the PR that 3.2.2 is sitting in
Same question as 2026-09-04 and 09-07, same answer: `$SNAKEMAKE_PROFILE` is a
new user-facing way to decide where a run executes, `--help` and `docs/`
changed with it, and a patch number would say "nothing new here" about the
thing the release exists to ship. The TUI's false-success fix on its own would
have been a patch.

**Pushing it today is free, and waiting is not.** Bioconda
[PR #68904](https://github.com/bioconda/bioconda-recipes/pull/68904) — "Update
comparem2 to 3.2.2" — has been open since 2026-09-07T21:25:46Z, so mergify's
`created-at<3 days ago` makes it eligible at 2026-09-10T21:25:46Z. The bot
pushes a newer tag onto that same `bump/comparem2` branch and rewrites the
title, so 3.3.0 inherits that clock and publishes at the same moment 3.2.2
would have. Tagging after 09-10 21:25Z instead would publish 3.2.2 and start a
fresh three days for the next number.

So **3.2.2 will never be published**, overwritten exactly as 3.1.0 and 3.2.0
were, and the channel goes 3.2.1 -> 3.3.0. Nothing is lost: 3.2.2's content is
in 3.3.0, and `run_exports` pins at `max_pin="x"` either way.

No recipe change: the diff added no dependency, no build-script change and no
test-section change, which are the only three things autobump cannot do.

### The docs ship a real report, and it keeps the section that contradicts itself
Asked for: a complete HTML report of an interesting genome set, in `docs/`, to
show new users what the product is. The set chosen is eight complete
*Streptococcus mitis* group genomes — six pneumococci across serotypes and both
pandemic MDR lineages, plus *S. pseudopneumoniae* IS7493 and *S. mitis* B6 — and
it is on the *Streptococcus mitis* group because that is where the 95% ANI
species boundary is genuinely contested, so three of the fourteen sections have
something to disagree about rather than merely reporting one species eight times.

**Including D39 *and* R6 was the load-bearing choice, and it is what caught the
defect.** R6 is D39's unencapsulated laboratory derivative; the run puts them at
100.00% ANI, 65 core SNPs, the same ST595 and branch lengths of 1.7e-5 and
4.2e-5 from their common ancestor. Twelve sections agree they are one strain. The biosynthesis panel calls **18 of 32
compounds *de novo* for D39 and 0 for R6**, flipping 20 of the 32 one way. That
is not a subtle disagreement and a pneumococcal reader would see it immediately.

**The set was chosen partly on a heuristic that then failed.** The reasoning was
that ~2,000-protein genomes sit below the ~3,200 where the *S. aureus* run of
2026-09-04 started hitting CarveMe's 600 s MILP ceiling. Measured here: **3 of 8
hit the ceiling anyway**, at gaps of 2.3–4.0% — *worse* than that run's 1.0–2.1%.
Proteome size does not predict it. Do not re-derive the heuristic; it is retired.

**And it is not a time limit away — measured, 18x.** `Spn_D39` re-solved at
`limits/time=10800` ran the full 3 h, returned `timelimit` at a 2.01% gap, and
came back with **1,135 reactions against 1,134** — one more, for eighteen times
the budget, still 444 short of R6's 1,579. `Smitis_B6` gained five. The sparse
model is therefore not a truncated version of the dense one; it is where the
search lands and stays.

**Nor is it about convergence.** `Spn_P1031` stopped on CarveMe's own gap
criterion (7.6e-4, under `limits/gap=0.001`), re-solved to a certified optimum at
the same objective 790, and *still* answers 18 de novo against 0 for four other
converged models. `Spn_R6` reproduces to the same objective 813.6 in 15.2 s, so
the problem is deterministic, not noisy. Two certified-optimal models, opposite
panel verdicts: this is the degeneracy STATUS.md already records for the
*S. aureus* set, now shown to reach the **verdicts** and not only the reaction
counts. That is a stronger statement than the one STATUS.md flagged as unmeasured
— the truncated model does not merely drop reactions, it *gains* de novo calls.

**Three options were considered for the showcase.** Omitting the two metabolism
sections was rejected: they are two of the three things `docs/index.md`
advertises, and a showcase that quietly drops the inconvenient section is worse
than one that explains it. Choosing a set without the near-identical pair was
rejected for the same reason and worse — it would have hidden the defect rather
than fixed it, and the pair is what makes any of this checkable. Fixing the
disclosure first (threading SCIP's status into the carveme and biosynthesis
sections) is the right change and is *not* done here, because `report.py` would
have to read a rule's log, which is not a declared output — that is a
`DESIGN.md` question, not a docs task.

So the report ships whole, and `docs/07 an example report.md` states the
contradiction with its numbers under a heading that says so. The page uses it to
teach the habit the repo already runs on: **put a known duplicate or
near-duplicate in your set and see which sections agree.** That is the standing
cross-check, handed to the user.

### GTDB-Tk's ANI-against-own-radius colouring earns itself on a three-species set
First time the section has rendered on genomes that disagree below genus. The
shared lineage collapses to `Bacteria › Bacillota › Bacilli › Lactobacillales ›
Streptococcaceae › Streptococcus` and one Species column carries
*S. pneumoniae* (6), *S. pseudopneumoniae* and *S. mitis_AR* — GTDB's split of
*S. mitis*, at 100.0% ANI because B6 is that cluster's own reference. All eight
green against a 95.0 radius at 97.53–100.0.

The point is what skani says about the same genomes: *S. pseudopneumoniae* sits
at **93.98–94.51%** ANI to the pneumococci, i.e. below the conventional 95%, while
GTDB-Tk still assigns it a species — because it is 97.53% to *its own* closest
reference. A report that coloured against a global 95% would have shown a
contradiction that is not there.

### The "0 de novo" models cannot take up nitrogen, and that was diagnosable without a biological assumption
Follow-up to the entry above, same day. The showcase left one thing open: which
side of the D39/R6 biosynthesis split was wrong. Answering it looked like it
needed a judgement about pneumococcal physiology, which is exactly the kind of
call that should not be guessed. It did not.

**The first hypothesis was wrong, and is kept because it was.** 2-oxoglutarate
looked like the hinge: no oxidative TCA cycle in a lactic acid bacterium means
no route from glucose, no glutamate, and no transaminated amino acid — one
missing link, twenty zeros. Measured, `akg` and `glu__L` are indeed unreachable
from M9 in **all eight** models. But that is what the two models *agree* on, so
it cannot be what separates them. A hypothesis that explains both sides of a
disagreement explains neither.

**What separates them is upstream of every enzyme:** four of the eight models
are missing a link in the three-step ammonium uptake chain, and M9's only
nitrogen source is ammonium. The partition against the panel is exact — the four
with `EX_nh4_e`, `NH4tex` and `NH4tpp` answer 17–18 de novo, the four without
answer 0. All eight carry `GLUDy`, `ASPTA` and `ALATA_L` with identical bounds.

So the zeros are defects, the 17–18 are the defensible numbers, and the working
route is gene-associated rather than gap-filled (`R_GLUDy` ← `G_AJOIJO_01205`,
with glutamate and 2-oxoglutarate cycling catalytically, which is why aspartate
is producible while glutamate is not).

**The generalisation, and the change it produced.** A draft metabolic model has
to be checked for whether it can take up the medium's carbon and nitrogen
sources *before* any phenotype is read off it. `biosynthesis.py` already
computed the ingredients — the media table's `present` column counts how many of
a medium's compounds the model has an exchange for — but it is a **count**, and
that is why this hid: R6 and D39 both read 17 of 20 for M9, differing only in
which three (`na1 nh4 ni2` against `mobd na1 ni2`), and 17 of 20 looks fine.

Now `SOURCES` names M9's one carbon and one nitrogen source,
`unreachable_sources()` probes each, the media TSV carries `missing` and
`unreachable` columns, and `_starved_note()` renders the sentence **directly
beneath the de novo counts** rather than as a footnote under the media table —
the qualification belongs next to the number a reader would quote.

**`unreachable` is a flux probe, not set membership on `probe.exchanges`, and
that is the load-bearing part.** The obvious implementation is the cheap one:
is the source compound in the exchange set? It would have cleared
`Spn_ATCC700669`, which reads the *highest* `present` count of the eight
(18 of 20), carries `EX_nh4_e` and `NH4tex`, and has no
periplasm-to-cytoplasm `NH4tpp` — so ammonium reaches the periplasm and stops.
It was also the model that reached a certified optimum with the most reactions
in the set. Three of the four broken models would have been caught by the cheap
check and the fourth, the most deceptive one, would not. A unit test pins the
behaviour with a stub probe rather than the implementation's text.

Only the M9 family gets the claim. LB carries nitrogen in every amino acid, so
no single compound is its source and a missing `nh4` there is not a defect;
asserting otherwise would have produced a false alarm on every rich medium.

The showcase report was re-rendered over the fix, so the shipped artefact
diagnoses its own limitation instead of relying on the docs page to do it. The
diagnostic scripts are in `upstream/` so the next session does not re-derive
them.

## 2026-09-09 — the logo goes in, and the theme's avatar slot has to be undone

The mark is three rows of contigs — genomes of differing length and
fragmentation — in the report's own palette (`#2b6cb0` accent, `#7fb3e8` tint).
`docs/assets/logo/make_logo.py` renders it: `comparem2-icon{,-dark,-badge}` and
`comparem2-logo{,-dark}`, each as SVG and PNG, plus the three concepts that
lost (synteny, dot plot, pangenome) which the script still draws.

**The rendered files are committed, not generated at build time.** The script
needs `cairosvg`, `fontTools` and Inter SemiBold, and `FONT` points at
`/tmp/inter/extras/otf/Inter-SemiBold.otf` — a path that will not survive a
reboot. Read the Docs installs `docs/requirements.txt`, which is mkdocs and
markdown-include and should stay that small. So the SVG/PNG are the artefact;
the script is the record of how they were made, and re-running it means
re-fetching Inter first.

**`docs/assets/extra.css` is load-bearing, not cosmetic.** The built-in
`readthedocs` theme styles its logo slot for a square avatar — theme.css pins
`.wy-side-nav-search img` to `45px` square with `border-radius:100%` on a
`#2980b9` backdrop. The wordmark is 331.7×56, a 5.9:1 lockup in near-black ink.
Dropped in unstyled it is squashed into a circle; on the blue header the `M2`
in `#7fb3e8` sits at about 2:1 contrast. Hence: un-crop the image, and make
that panel white so dark ink has something to sit on. Deleting the file does
not return the sidebar to the theme default — it returns it to a cropped
wordmark.

**The favicon is `docs/img/favicon.ico` because that is the only hook left.**
`site_favicon` was removed from mkdocs' config schema (it is absent from
`config/defaults.py` in 1.6, so `--strict` would reject it), but both built-in
themes still emit `<link rel="shortcut icon" href="img/favicon.ico">` as the
fallback, and a file under `docs/` overrides the theme's copy. Verified against
the build, not assumed: `site/img/favicon.ico` is byte-identical to ours. The
ICO carries 16–256 px from the badge variant — the blue-square one, because at
16 px a transparent mark on browser chrome is three faint bars and the badge is
still recognisably the mark.

**The homepage leads with the wordmark instead of `# CompareM2`,** which is why
`nav` now says `Home: "index.md"` — the label used to come from that heading.
Two copies of the same lockup on one screen (sidebar and hero) was the
alternative and it read worse.

**The README uses `<picture>` with `prefers-color-scheme`** and relative paths.
Both were checked rather than assumed: GitHub's `POST /markdown` keeps
`<picture>`/`<source>` through sanitisation and wraps them in
`themed-picture`, so the dark variant does get used. Relative paths are right
because nothing renders this README off GitHub — `readme = "README.md"` in
`pyproject.toml` feeds a PyPI long description that is never published, since
the bioconda recipe builds from a GitHub tag tarball.

`exclude_docs` now keeps `generate.py`, the logo script and the 315 kB concept
sheet out of the built site. They are inputs to the documentation, not pages of
it.

**The report gets the mark too, and it is drawn rather than embedded.** The
heading carries it as inline SVG whose rects are filled from `var(--accent)`
and `var(--tint)` — a new custom property, `#7fb3e8` light and `#3f80c4` dark,
which is what `make_logo.py` already used for the two icon variants. Painting
from properties is the whole point: an `<img>` or a `background-image` cannot
see them, so the mark would sit in light-mode blue on a `#161616` page.
Measured on the two schemes: accent 5.42:1 on white and 8.20:1 on `#161616`,
tint 2.21:1 and 4.39:1. The tint is the quiet half of a decorative mark next to
a text heading, not information, so the low light-mode figure is the design and
not a defect.

`h1 .mark` has to set `display:inline-block`. Every other SVG in the report is
a full-width figure, so the sheet sets `svg { display:block }`, and the mark
first rendered on its own line above the title — 29.4 px in the right place,
one line too high. Cheapest thing to get wrong twice.

**The favicon is a `data:` URI holding the badge variant with literal
colours**, which cost one narrowing of `test_report_is_self_contained`: it
asserted `"<link" not in body`, and the rule it means is *nothing the browser
fetches*, so it now allows a `<link>` whose href is a data URI and fails on any
other. Literal colours because a data URI is a separate document and cannot see
the page's properties — the one place the mark cannot follow the reader's
scheme. The badge, not the plain mark, because at 16 px a transparent mark on
browser chrome is three faint bars of unspecified colour; blown up 10x, the
badge still reads. Safari may ignore an SVG favicon entirely and show its own
default; that is untested here.

The geometry is duplicated between `report.py` and `make_logo.py` rather than
shared. `docs/` is not in the package — `plasmids.zip` is the only non-Python
file that ships — so importing or reading the artwork would either break that
or add package data for five rectangles.

`docs/assets/example-report.html` was **patched, not re-rendered**: the run it
came from is on GenomeDK. The substitution is provably what the current
renderer emits, because the example's inlined `<style>` was byte-identical to
`report.py`'s `CSS` constant beforehand — so replacing that block, adding the
favicon `<link>` and putting the mark in the one `<h1>` reproduces today's
output for that run exactly, and the diff is five hunks all inside the head.
156,180 → 158,052 bytes. If a future change to `CSS` lands while the example
still needs updating, the same three substitutions are the recipe; a real
re-render needs the cluster.

### The prose drops the em dash, because the tells were mechanical and the content was not
Carl pointed at
[WP:AISIGNS](https://en.wikipedia.org/wiki/Wikipedia:Signs_of_AI_writing) and
asked for documentation that does not read as machine-written. Measured first:
309 em dashes across README and `docs/`, 17.4 per 1,000 words, against the 1–3
that ordinary technical prose runs. Also 52 rhetorical `rather than`
antitheses and 36 `**Bold lead-in.** Explanation` paragraphs. Now 28 and 1.5.

What the essay actually weights most was already absent, and that is the part
worth recording: no puffery vocabulary, no `serves as` for `is`, no "not only X
but Y", no vague attribution, no invented citations. Checked by grep after the
rewrite as well as before. So this was a surface pass and the content did not
move — verified by diffing every numeric token in `guidance.py` against the
previous commit (one `82%` gained a comma where a dash had been; nothing else)
and by comparing all 25 verbatim paper quotations byte-for-byte.

**The single biggest source was one line of `docs/generate.py`**, not the
prose. `parts.append(f"  - *{label}* — {text}")` emitted roughly 65 of page
30's em dashes on its own; it is now a colon. Worth knowing before anyone
hand-edits that page: it is generated, and the separator lives in the
generator. The report never had this problem, because `report.py` renders the
same pairs as a real `<dl>`/`<dt>`/`<dd>`.

Kept deliberately: em dashes inside verbatim UI output (`databases to download:
none — all present`, the quit dialog, the spinner's idle line), inside `| — |`
table cells meaning "not applicable", and in the FastTree 2 paper's actual
title. Those are quotations and data, not voice. Also kept: 21 `rather than` on
page 30, because it is a plain comparative that the essay does not flag, and in
most of them it carries the claim — "a statement about the model rather than
the organism" does not survive being reworded to hit a count.

Three stale test counts were fixed in passing. `README.md` said 276,
`docs/10 installation.md` 269, `CLAUDE.md` 275; CI and a local run both say
**274 passed, 2 skipped**.

## 2026-09-09 — the sidebar wordmark shrinks, and the white panel is reversed

Same day as the logo went in, and it undoes that entry's last move. Making the
sidebar header white so dark ink had something to sit on solved the contrast
problem and created a worse one: on a narrow screen the header sits directly
under `.wy-nav-top`, the theme's blue mobile bar, and the two backgrounds met
at a hard edge. Carl saw it on a phone. The panel was also carrying a 250 px
wordmark, so the same lockup appeared twice on the homepage at nearly the same
size.

The header keeps blue and the wordmark changes instead. `make_logo.py` now
renders a third lockup, `comparem2-logo-badge`, from the `badge` palette it
already had for the icon — white ink, `#bcd8f5` for the `M2` and the middle
contig row. That palette needed one new key, `wordmark_accent`, because
`badge` sets `accent` to white for the mark and the `M2` would otherwise
disappear into `Compare`. **The lockup leaves `bg` unused**: the blue comes
from whatever it sits on, so the header and the mark cannot drift apart.

Measured on the header blue, `#2b6cb0`: white 5.42:1, `#bcd8f5` 3.69:1 — above
3:1, and the wordmark is large text (cap height 15 px at the rendered size,
SemiBold). The reason the existing dark-background lockup could not be reused
is the same table: its `#7fb3e8` `M2` reads **2.46:1** there, and its `#3f80c4`
contig row is nearly invisible.

Both blue bars are now `#2b6cb0` rather than the theme's `#2980b9`, so the docs
chrome and the report agree on the accent, and white-on-blue goes from 4.30:1
to 5.42:1.

Size: 27 px tall, flush left, aligned to the search field's left edge (a −4px
margin cancels the link's own padding). Header height 101 px against 128 px
before. The homepage hero stays at 420 px — it is the page's title, and the
duplication was two *large* copies, not two copies.

**The theme's rule for that slot is `.wy-side-nav-search > a img.logo`.** A
plainer `.wy-side-nav-search img.logo` loses on specificity no matter how late
`extra.css` loads, and the first attempt at this change silently kept the
theme's `height:auto` — the header looked untouched. If a size here stops
taking effect, that is where to look. The `:hover` override went the other way:
it is *removed*, so the theme's white-10% highlight comes back now that there
is a blue background for it to lighten.

Verified in a real build — `mkdocs build` into a temp dir, screenshotted at
1280 and 390 px wide with the mobile nav both closed and open — rather than
read off the CSS. mkdocs is not a project dependency, so that was a throwaway
venv from `docs/requirements.txt`.
