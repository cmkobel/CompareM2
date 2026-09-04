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
