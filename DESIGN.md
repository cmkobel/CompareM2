# CompareM2 v3 — the design

What v3 is and why it is shaped this way. No dates and no history here: the
dated record of how it got this way is [DECISIONS.md](DECISIONS.md), and what is
currently true of a real run is [STATUS.md](STATUS.md).

Only one thing survives from v2: **easy to install, easy to run, easy to
interpret** — a straightforward view into a set of assemblies.

## Wide, not deep

v3 is triage across many assemblies, not a deep dive into each. Deep functional
and metabolic analysis is where DRAM2 and nf-core/funcscan are already strong; a
fast wide sweep is not well served, and it is what v2's benchmark result —
near-linear scaling with input count — actually supports.

Every other decision follows from that one. Thirteen tools instead of 30+, one
conda environment instead of 25, one runtime instead of two.

## The shape

**Tools are declared once; the workflow is derived.** There is no hand-written
Snakefile. `catalogue.py` holds thirteen `Tool` specs and `snakefile.py`
generates a Snakefile from them, so adding a tool cannot mean editing a
workflow.

```
src/comparem2/
  tools.py      the contract: Tool, Database, Context, Registry, Scope
  catalogue.py  the thirteen specs — command lines and outputs
  guidance.py   what each tool does and how to read it, for the report
  snakefile.py  generates the Snakefile and per-tool env files
  cli.py        arguments, input canonicalisation, hands off to Snakemake
  runner.py     drives Snakemake's API, emits structured events
  tui.py        Textual interface over those events
  report.py     renders the HTML report
```

The trick that makes generation work: `Context` renders every path from a sample
name, so building one with `sample="{sample}"` and `threads="{threads}"` yields
Snakemake wildcards from exactly the same code that yields real paths at report
time.

**Snakemake is kept, as an executor rather than a framework.** The DAG was never
the hard part — *cluster execution* is: sbatch generation, dependency chains,
`squeue` polling, partial-failure recovery, retries. Snakemake's SLURM executor
plugin already solves it, and it is tens of MB against 141 GB of databases.

**Python only.** CLI, TUI and report share one runtime. R, pandoc, rmarkdown and
tidyverse are gone — they were the heaviest non-database dependency in v2.

**A tool declares its outputs** rather than writing wherever it likes. That
single choice buys resumability, progress reporting, and the report's knowledge
of what to render.

## Install weight is a first-class constraint

Every database declares its size in bytes, **measured** from `content-length`,
and `cm2` shows what it is about to fetch before fetching it. `None` means
unmeasured and is reported as such — a total that silently omits an unknown is
worse than no total.

The environment is slim; the data is not, and GTDB-Tk is the whole reason. At
141.4 GB it is 91% of the install. Making it opt-in was offered twice and
declined twice; it stays in the default path. If "easy to install" ever has to
be defended, that is the line item.

**One environment, plus CheckM2 on its own.** CheckM2 pins DIAMOND 2.1.x while
current Bakta needs 2.2.x, so they cannot co-solve. CheckM2 is the single tool
marked `isolated=True`. Everything else shares one environment.

**Databases fetch themselves, as workflow rules.** Each `Database` declares a
`fetch` and a `ready` path, and a `download_<name>` rule is generated whose
output is that path. Every tool takes its database's `ready` path as an input.
Downloads therefore inherit what Snakemake is already good at: a present
database is skipped, a half-finished one is redone, fetching overlaps unrelated
work, and ordering is the DAG's problem. A separate `--download` phase would
have reimplemented all of it.

`ready` is a real file the tool needs wherever possible, not a stamp, because a
stamp can outlive the data it claims.

## The report is the product

Once the tools are commodity, the report is what is left. Two rules:

**Every tool gets a section.** A generic table renderer means a tool cannot
silently produce nothing — which is how v2 ended up running antiSMASH,
InterProScan, IQ-TREE, FastTree and TreeCluster without displaying any of them.
A specific renderer is an improvement on the fallback, never a prerequisite.

**Sections appear only when their outputs exist**, so a partial run still
produces a readable document. Output is one self-contained HTML file with no
external assets, so it survives being emailed or copied off a cluster.

**The report explains itself.** Every tool carries a `Guidance` value: what
question it answers, how it works, what each rendered column means, and what the
result cannot tell you — as a collapsed block under each heading, plus a
"Methods and citations" list covering exactly the tools that ran.

Guidance lives outside `catalogue.py` because a tool's spec there is ~20 lines
of command line, and burying it under prose written for a different audience
hides the thing developers edit. The single-source-of-truth invariant is kept by
a test instead: `test_every_tool_has_guidance` fails when a tool is added
without an explanation.

**Two docs pages are generated** from the specs — the analyses reference and the
citation list — and CI fails if they are stale. The analyses page renders each
tool's *actual* command line, so it cannot describe something that does not run.

## Rules that must not be quietly undone

Each of these exists because something broke, or because undoing it would break
something already published. The post-mortems are in
[DECISIONS.md](DECISIONS.md#what-went-wrong); the short version:

- **The name stays `CompareM2`, whatever the version number says.** The 2 is not
  a version, so `CompareM2 v3.0.0` is correct even though it reads oddly, and
  renaming to match the major version is the obvious-looking move that must not
  be made. The citation trail is fixed and cannot be renamed with the code: the
  paper is doi:10.1093/bioinformatics/btaf517, it names CompareM2, and it sends
  readers to `comparem2.readthedocs.io`. The bioconda package v2 shipped under
  is the same string. A rename orphans all three and makes a published tool look
  like a different one — and the repository has already absorbed one rename
  (`comparem2` → `CompareM2`), which still emits redirect notices.
- **Pin a minimum version for every tool.** An unconstrained bioconda spec lets
  the solver reach back years to satisfy some other package's constraint. Two
  tools resolved to builds that installed cleanly and crashed on first use.
  **`pixi install` succeeding says nothing about whether the pipeline works.**
- **Commands are argument lists, never shell strings.** A tool that writes to
  stdout declares `stdout_to_output=True`; the redirect is added by whatever
  runs it. The same discipline applies to database `fetch` steps.
- **`isolated=True` is an exception that must carry its reason in the spec.**
  v2 reached 25 environments by making it the default.
- **A tool's database location must be reachable.** Some tools take it only
  through the environment, which is what `Tool.env` is for — without it,
  `--databases` was silently ignored for the largest database in the pipeline.
- **The database directory defaults to a shared, home-relative location**
  (`~/.comparem2/databases`, or `$COMPAREM2_DATABASES`), never one relative to
  the cwd or the output directory. Databases outlive any one run's results, and
  a per-run default silently buys a second copy of 143 GB.
- **Unit tests are the primary instrument.** The codebase is a generator, and a
  wrong wildcard yields a Snakefile that parses cleanly and builds the wrong
  DAG. An end-to-end run catches that slowly, if at all.
- **Verification tracks execution, never installation.** A drafted command line
  is worthless until it has been run on real genomes.
- **Numbers in `guidance.py` are quoted from the paper and checked against it.**
  Confabulation dressed as insight is the failure mode that reads best, so it
  gets a deterministic gate rather than a second opinion.

## Open questions

- **Taxonomy.** The single largest install cost. GTDB-Tk is authoritative and
  141.4 GB; skani or mash against GTDB are fast and small but approximate.
  Candidates not yet measured. Whatever replaces or precedes GTDB-Tk **must take
  assemblies** — that was the mistake sylph embodied.
- **GTDB-Tk's summary output.** It writes `bac120` and `ar53` summaries
  separately; the declared `gtdbtk.summary.tsv` needs a concatenation step that
  does not exist yet. GTDB-Tk has also never been run.
- **CarveMe should probably be opt-in.** It works, but at roughly nine minutes
  per genome it is by far the slowest per-genome tool, which matters for a wide
  view over hundreds of assemblies.
- **AMRFinder's database cannot live under `--databases`.** `amrfinder -u`
  rejects `-d`, so its data goes into the conda prefix and all that can be
  recorded is a stamp — which does not survive the environment being rebuilt.
- **FastTree runs single-threaded.** Parallelism needs the `FastTreeMP` build
  plus `OMP_NUM_THREADS`; `Tool.env` now makes that possible but it is not done.
- **Report sections** still to write from scratch: v2 never displayed fasttree or
  treecluster, so those are new rather than ports.
