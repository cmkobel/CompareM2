# CompareM2 v3 — the design

What v3 is and why it is shaped this way. No dates and no history here: the
dated record of how it got this way is [DECISIONS.md](DECISIONS.md), and what is
currently true of a real run is [STATUS.md](STATUS.md).

Only one thing survives from v2: **easy to install, easy to run, easy to
interpret** — a straightforward view into a set of assemblies.

## Wide, not deep

v3 is triage across many assemblies, not a deep dive into each. Deep functional
and metabolic analysis is where DRAM and nf-core/funcscan are already strong; a
fast wide sweep is not well served, and it is what v2's benchmark result —
near-linear scaling with input count — actually supports.

Every other decision follows from that one. Fourteen tools instead of 30+, two
conda environments instead of 25, one runtime instead of two.

## The shape

**Tools are declared once; the workflow is derived.** There is no hand-written
Snakefile. `catalogue.py` holds fourteen `Tool` specs and `snakefile.py`
generates a Snakefile from them, so adding a tool cannot mean editing a
workflow.

```
src/comparem2/
  tools.py      the contract: Tool, Database, Context, Registry, Scope
  catalogue.py  the fourteen specs — command lines and outputs
  guidance.py   what each tool does and how to read it, for the report
  snakefile.py  generates the Snakefile and the two env files
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
plugin already solves it, and it is tens of MB against 62.5 GB of databases.

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
60.8 GB it is 97% of the 62.5 GB the pipeline measures. Making it opt-in was
offered twice and declined twice; it stays in the default path. If "easy to
install" ever has to be defended, that is the line item.

It was 141.4 GB until 2026-09-02, and the reduction was not a choice: GTDB-Tk
2.7 requires r232 and refuses r226, and r232 dropped FastANI's reference
genomes in favour of skani sketches. **A tool's database version is part of its
pin.** The catalogue had been carrying an r226 URL that the installed tool
would have rejected — after the download.

**Two environments, and CheckM2 is why there are two.** It pins DIAMOND 2.1.x
while current Bakta needs 2.2.x, so they cannot co-solve. Everything else shares
`main`. Snakemake deploys both; nothing expects a tool on PATH.

**Databases fetch themselves, as workflow rules.** Each `Database` declares a
`fetch` and a `ready` path, and a `download_<name>` rule is generated whose
output is that path. Every tool takes its database's `ready` path as an input.
Downloads therefore inherit what Snakemake is already good at: a present
database is skipped, a half-finished one is redone, fetching overlaps unrelated
work, and ordering is the DAG's problem. A separate `--download` phase would
have reimplemented all of it.

`ready` is a real file the tool needs wherever possible, not a stamp, because a
stamp can outlive the data it claims.

## Steps around a command

Eleven of the fourteen tools are one command: arguments in, declared files out.
GTDB-Tk is not, and rather than let it become a hand-written rule it gets two
fields that every tool could use and only it does.

- **`Tool.files`** — input files the tool needs that are nobody's output, as a
  mapping of path to content rendered from the `Context`. `prepare()` writes
  them and the rule declares them as inputs. GTDB-Tk needs a two-column
  `--batchfile` because canonicalisation puts each genome in its own directory,
  so `--genome_dir` cannot be used.
- **`Tool.post`** — argument lists run after the command, for turning what a
  tool writes into what its spec declared. GTDB-Tk writes `bac120` and `ar53`
  summaries separately and the pipeline declares one table.

CarveMe needs the third shape: a step *in front of* the command rather than
around it. `carve_scip.py` is the whole of it — it turns off one SCIP presolver
and hands `carve` an input path inside CarveMe's own output directory, then
calls CarveMe in-process. Both halves are things only the calling process can
do: the solver parameter because neither CarveMe nor ReFramed exposes one, the
input path because `carve` derives its DIAMOND output from it and would
otherwise overwrite Bakta's feature table.

`biosynthesis` is the fourth shape and the only tool here that *is* ours:
`biosynthesis.py` is the program, not a wrapper around one. It reads CarveMe's
model with ReFramed and answers, per compound, whether the network has a route
to it. That it is a tool rather than a report section is deliberate — it has a
declared output, so it is resumable, it appears in `--until`, and the report
reads a file instead of running a solver.

Both run under a **bare `python`**, which is the exact opposite of the rule
below, and for the mirror-image reason: they import `carveme` and `reframed`, so
they need the interpreter of the environment the *tool* is in. `steps.py` is our
code and needs ours. All three import nothing from their own package, because
none of them run in an environment that has it — and
`biosynthesis.py` additionally keeps its `reframed` imports inside the functions
that solve, so `report.py` can read the compound panel from it in a process
that has no solver.

Three properties matter more than the fields:

**They are declared, not hidden in a shell string.** The batchfile is a rule
input, so the DAG knows about it. Before this, the rule named a batchfile that
nothing created and declared an output the tool never writes — both described
in comments as though they existed, in a rule that had never been executed.

**Written only when changed.** A declared file is a rule input, so rewriting an
identical one moves its mtime and re-runs the tool that reads it. For GTDB-Tk
that is hours of work triggered by a four-line file.

**The step runs through an absolute `sys.executable`.** Under a conda
deployment the rule's environment holds the tool, not CompareM2, so a bare
`python -m comparem2.steps` would find the wrong interpreter.

`Tool.post` is not a licence for shell work: a tool needing several steps is a
tool whose spec is lying about what it does, and a test asserts that GTDB-Tk is
the only tool using either field.

## One deployment model, and six environments

**Snakemake installs the tools. That is the only way a tool arrives.** There is
no flag, no fallback, and no mode in which the pipeline expects a tool on PATH.
`pixi` is how a developer gets *the pipeline* — python, Snakemake, textual,
pytest, conda — and it is deliberately not how anyone gets the tools.

This replaced a two-model design on 2026-09-03, in which pixi solved one
environment holding thirteen tools and `--use-conda` switched to per-rule
deployment. Two models meant the tool set was pinned in two files that could
drift, and they had: `gtdbtk>=2.7` was pinned in `catalogue.py` and left
unconstrained in `pixi.toml`. It also meant a conda-installed user's first run
ended in `not on PATH: seqkit ...` and an incantation to remember. Both are
gone; see DECISIONS.md.

**Eighteen rules, six environments, grouped by dependency ecosystem.**
`catalogue.py` declares them, and a test enforces the list:

- `basic` — seqkit, skani, snp-dists, fasttree, treecluster, plus curl and tar.
  Small standalone binaries with no deep dependency tree, which is why the two
  database fetches that are a plain URL run here too.
- `perl` — mashtree, mlst, panaroo. The fault line; see below.
- `annotation` — bakta, amrfinder, and both of their database fetches.
- `gtdbtk` — one tool.
- `carveme` — carveme and biosynthesis, which imports `carveme` and `reframed`.
- `checkm2` — one tool, for one reason. CheckM2 pins DIAMOND 2.1.x and a
  current Bakta needs 2.2.x, so an environment holding all fourteen does not
  exist. Verified 2026-09-01: `bakta>=1.10` solves with all thirteen others and
  fails only when checkm2 joins.

**This was two until 2026-09-07, when the thirteen-tool `main` stopped solving
at all.** Not on one machine and not from new configuration: the identical spec
built a working 6.0 GB environment on thylakoid on 09-03 and, re-solved there on
09-07, failed after 6 min 46 s — and after 4 min 07 s on GenomeDK. The Perl
stack went bad upstream and took the other ten tools with it. Every group above
solves on its own in under 30 s: `basic` 12 s, `annotation` 14 s, `gtdbtk` 17 s,
`carveme` 14 s, `checkm2` 26 s, `perl` 28 s, measured the same day on GenomeDK.

**A co-solve is a shared fate.** Thirteen tools in one environment means
thirteen sets of transitive constraints that must hold at the same moment, so
any one ecosystem going bad takes the rest down — and it did. Grouping by
ecosystem is what keeps one bad month in bioconda from making the whole pipeline
uninstallable, and it is what two was trading away for a disk saving.

**An environment per *tool* is still wrong.** It is what content addressing
makes easy and it is v2's 25 environments in another form: fourteen solves cost
fourteen copies of DIAMOND, python and numpy, and the user pays on first run.
Six is the number of distinct ecosystems, not a compromise between one and
fourteen — and per-rule deployment remains the only distribution model that can
work at all, because a recipe listing the tools as run dependencies would still
have to drop CheckM2 or Bakta.

Six consequences worth knowing:

- **Environments are addressed by content, not by rule.** Snakemake deploys to
  `md5(realpath(conda_prefix) + env file content)`, so `render_envs` writes one
  file per *environment* and eighteen rules point at six of them.
- **The prefix is in the identity, and the workdir is not.** One `basic.yaml`
  deployed to three prefixes gets three directories, so `--conda-prefix` is the
  one argument a later run has to match; but `--output` may differ freely,
  which is what makes `--setup` able to build in a temp directory that it then
  deletes. Both halves measured 2026-09-03.
- **`--setup` exists because the build happens before the first job**, not
  lazily per rule — so a first run is silent for as long as the solves take,
  which on a cluster is when someone kills it. `cm2 --setup` asks for the same
  work on purpose, and needs neither assemblies nor databases: a stub FASTA
  closes the DAG, and every database path in it is some rule's output.
  **It cannot be done at install time**, which was asked and answered: the
  prefix is a *runtime* choice, so a package that deployed to the default
  location would have built the wrong directory for anyone using
  `$COMPAREM2_CONDA_PREFIX` — on a cluster, everyone.
- **Splitting narrows the blast radius; it does not remove it, and the pins are
  load-bearing.** Grouping by ecosystem means a bad Perl month costs three tools
  rather than thirteen, but every group is still a co-solve. An unconstrained
  spec inside one of them reaches back years to satisfy some other package's
  constraint, so **every tool carries a minimum version**, floored at the build
  verified on linux-64, and a test enforces it. The versions every verification
  run on thylakoid ran on: seqkit 2.13.0, bakta 1.12.1, panaroo 1.8.0,
  gtdbtk 2.7.2, DIAMOND 2.2.5.
- **AMRFinder depends on the sharing.** `amrfinder -u` writes into
  `$CONDA_PREFIX`, so its download rule and its analysis rules must land in the
  same deployed environment. They do — both name `annotation`, as `main` before
  it. **Executed 2026-09-03**, when they were still separate single-tool
  environments sharing a spec string: five rules, one directory, the database
  inside it, and AMRFinder's own log naming that directory as both its software
  and its database path. This is a constraint on any future regrouping: bakta
  and amrfinder may move, but not apart. See STATUS.md.
- **A database fetch is a rule and needs an environment too.** Two of the four
  run a tool binary rather than curl (`bakta_db download`, `amrfinder -u`, both
  in `annotation`) and the other two need only curl and tar (`basic`), so
  `Database` declares `conda` and `environment` exactly as `Tool` does. Under a
  cluster profile they are additionally `localrules` — a compute node with no
  outbound network cannot fetch 60.8 GB of GTDB.

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
- **A solver is part of the pinned surface too, and its output is not just
  slower or faster — it is different.** CarveMe's MILP took 601 s and returned a
  model missing 253 of its annotated reactions with the SCIP conda-forge ships,
  and 10 s with the one the PyPI wheel ships, on a byte-identical problem. The
  difference is PaPILO, which conda-forge links. Whether that is a presolve bug
  or a consequence of CarveMe's bigM=1e3-against-eps=1e-3 scaling is unsettled —
  what is measured is that this presolver costs 12–60x wall time and a quarter
  of the reactions the annotation supports. `carve_scip.py` turns it off;
  removing that line costs nine minutes a genome *and* the model. Re-measure
  before touching it — the numbers are in the wrapper's docstring and in
  [STATUS.md](STATUS.md).
- **A database version is part of the tool's pin.** Both directions are runtime
  failures that no solve catches: GTDB-Tk 2.7 accepts only r232, Bakta 1.12
  only db 6.x. Moving either one alone is the bug. Which is why the URL, the
  size and the `>=` pin sit in the same spec.
- **Commands are argument lists, never shell strings.** A tool that writes to
  stdout declares `stdout_to_output=True`; the redirect is added by whatever
  runs it. The same discipline applies to database `fetch` steps.
- **"Has this already run" is answered from the declared outputs, by one
  function.** `tools.completion()`, counted per unit of work — per genome for a
  genome-scope tool. The TUI's opening state, `any_outputs_exist` and the
  report's section-by-section rendering are the same question at three
  thresholds (`done`, `complete > 0`, `started > 0`), and they were three
  hand-written spellings before. Declared outputs are also what Snakemake
  decides resumability on, which is what makes the TUI's table a statement
  about what a re-run will skip rather than a guess. Results found on disk are
  a *different state* from results this session produced.
- **Every environment must carry its reason in `catalogue.py`.** There are six,
  each one an ecosystem: the DIAMOND conflict for `checkm2`, the 2026-09-07 Perl
  breakage for the rest, both written above the specs with the measurement that
  justified them. A test enforces the list. v2 reached 25 environments by making
  isolation the default rather than the exception, and an environment per *tool*
  is that mistake in a cheaper disguise — six is the count of ecosystems, and a
  seventh needs a reason of the same kind, not a preference.
- **A tool's database location must be reachable.** Some tools take it only
  through the environment, which is what `Tool.env` is for — without it,
  `--databases` was silently ignored for the largest database in the pipeline.
- **The database directory defaults to a shared, home-relative location**
  (`~/.comparem2/databases`, or `$COMPAREM2_DATABASES`), never one relative to
  the cwd or the output directory. Databases outlive any one run's results, and
  a per-run default silently buys a second copy of 143 GB.
- **The conda prefix defaults the same way** (`~/.comparem2/envs`, or
  `$COMPAREM2_CONDA_PREFIX`), for a sharper version of the same reason:
  Snakemake includes the prefix's realpath in each environment's hash, so
  moving it re-solves all six environments *and* re-fetches AMRFinder's
  database, which lives inside `annotation`.
- **Cluster submission is Snakemake's, and so is the variable that turns it
  on.** `--profile` defaults from `$SNAKEMAKE_PROFILE`, declared `env_var=` on
  Snakemake's own `--profile`, and there is deliberately no
  `$COMPAREM2_PROFILE` — v2 had one, and a second name for someone else's
  setting needs a precedence rule and buys nothing. Two consequences read like
  noise and are not: every Snakemake launched from `cli.py` names its profile
  explicitly, `none` included, because otherwise an exported variable submits a
  run this program has just reported as local; and the `--cores` gate keys on
  the *effective* profile, not on `args.profile`, which is what it did for as
  long as the variable went unread here.
- **The bioconda package must not grow tool dependencies.** It ships the
  pipeline; Snakemake deploys the tools. Adding them to the recipe would
  require dropping CheckM2 or Bakta — see *One deployment model*. The same
  rule applies to `pixi.toml`, which listed them until 2026-09-03 and drifted
  from `catalogue.py` while it did.
- **Snakemake is invoked as `sys.executable -m snakemake`**, never as a bare
  `snakemake`. It is this package's own dependency, so the correct one is the
  one beside the running interpreter; by name it was simply not found when a
  packaged `comparem2` was invoked by absolute path.
- **The report's own artwork is drawn in Python, never linked or embedded as a
  file.** The mark in the heading is inline SVG painted from `var(--accent)`
  and `var(--tint)`, so it inverts with the reader's colour scheme — an `<img>`
  or a `background-image` cannot see a custom property and would sit in
  light-mode blue on a dark page. The favicon is the exception that shows the
  reason: a `data:` URI is a separate document, does not see the page's
  properties, and so carries literal colours and cannot follow the scheme.
  Shipping either as an asset would break the other invariant too —
  `plasmids.zip` is the only non-Python file in the package, and `docs/` is not
  in it, so the geometry is duplicated in `report.py` on purpose.
  `test_report_is_self_contained` permits a `<link>` only with a `data:` href.
- **The interface indicates state with attributes, not with tints.** Every
  cursor and highlight is `text-style: reverse`, and the tempting cleanup — put
  the theme's own colour back — is the bug it was written for. tmux ships
  `default-terminal screen`, so eight colours is the ordinary case over SSH,
  and Textual's blurred cursor is `$primary` at 30% alpha: it blends to
  `#153854`, downgrades to ANSI 8 against a surface on ANSI 0, and disappears.
  The command palette showed it first because its list can never take focus.
  No colour can be made safe here — the two ANSI themes set `surface` to
  `ansi_default`, the user's own terminal background — while `reverse` is SGR 7
  and inverts whatever is already there. Verified at the byte level: a
  `standard`-system console emits `\x1b[7;36;40m`.
- **Stopping a run is done by this code, not by Snakemake.** Its scheduler
  reaches `executor.cancel()` only from a `KeyboardInterrupt` inside its own
  loop, and installs the SIGTERM handler that would get it there inside a
  `try/except ValueError` that skips silently off the main thread — which is
  where `runner.run()` puts it. So `cancel.py` cancels a queue with
  `scancel --name <run_uuid>`, the name the SLURM plugin submits under for
  exactly this purpose, and signals the process tree for a local run, because
  Snakemake's local `cancel()` is `pool.shutdown()` and waits rather than
  killing. Both halves run under a profile: the `download_*` rules are
  `localrules`, so a 60.8 GB fetch is a child process and not a job.
- **Unit tests are the primary instrument.** The codebase is a generator, and a
  wrong wildcard yields a Snakefile that parses cleanly and builds the wrong
  DAG. An end-to-end run catches that slowly, if at all.
- **Verification tracks execution, never installation.** A drafted command line
  is worthless until it has been run on real genomes.
- **Numbers in `guidance.py` are quoted from the paper and checked against it**,
  or, where a number is a measurement of this pipeline rather than a published
  result, the sentence says so. Confabulation dressed as insight is the failure
  mode that reads best, so it gets a deterministic gate rather than a second
  opinion.
- **A compound on the biosynthesis panel is probed in the form its pathway ends
  at, and no two members may come from one nutrient family.** Free thiamine,
  folate and lipoate are salvage substrates rather than biosynthetic products:
  probing `thm` called *E. coli* a thiamine auxotroph, probing `fol` called
  *S. aureus* unable to make the compound sulfonamides work by blocking, and
  probing `lipoate` returns "cannot make" for every organism there is. Two
  members of one family — `nac` and `nad`, say — rescue each other in the second
  test and the pair then reports a kinase rather than a pathway.
- **A panel target that is the same answer for every genome is not a target.**
  `btn` and `q8` were on the panel until 2026-09-10 and were `none` or `absent`
  in 12 of 12 CarveMe drafts across 8 species, `de novo` only in a curated
  model — including drafts of four organisms that make biotin and two that use
  ubiquinone-8. Both routes are in the reaction database and carving does not
  keep them, so the columns were two false dependencies a genome rather than a
  signal. All three rules are enforced by tests, and the calibration is now
  **29 of 30 for the curated `iML1515` and for CarveMe drafts of five
  prototrophs alike**, adenosylcobalamin the only miss in any of them.
  `upstream/panel_calibration.py` re-measures it from proteomes CarveMe bundles.

## Open questions

- **Taxonomy.** The single largest install cost. GTDB-Tk is authoritative and
  60.8 GB; skani or mash against GTDB are fast and small but approximate.
  Candidates not yet measured. Whatever replaces or precedes GTDB-Tk **must take
  assemblies** — that was the mistake sylph embodied. Less pressing than it was:
  r232 more than halved the download.
- ~~**GTDB-Tk's summary output.**~~ Closed 2026-09-02: `Tool.post` runs
  `comparem2.steps merge-tsv` after the command, and `Tool.files` writes the
  batchfile the rule had only claimed to write. See *Steps around a command*.
- ~~**CarveMe should probably be opt-in.**~~ Closed 2026-09-02: the nine minutes
  were a solver defect, not the method. 62 s a genome measured end to end
  through `carve_scip.py`, on one core, against Bakta's 77 s on eight for the
  same genome. Whether it should be opt-in anyway is a smaller question than it
  was, and no longer a performance one.
- **AMRFinder's database cannot live under `--databases`.** `amrfinder -u`
  rejects `-d`, so its data goes into the conda prefix and all that can be
  recorded is a stamp — which does not survive the environment being rebuilt.
- ~~**FastTree runs single-threaded.**~~ Closed 2026-09-03: `FastTreeMP` and
  `Tool.env` would make the switch reachable, and it was measured to be worth
  nothing — **including at hundreds of genomes**, which was the obvious
  objection and the reason the measurement went to 400 taxa. At 4 taxa MP is
  2.3% slower; at 7 every thread count lands inside the plain binary's own
  run-to-run spread. The speedup from 8 threads then does not grow with taxon
  count: **1.13x at 25, 1.15x at 100, 1.09x at 400**, for 40–65% more CPU and 8
  cores Snakemake would reserve. `threads=1` is a measurement now, not an
  assumption.
- **The support phase is most of FastTree's runtime at the shape this pipeline
  actually runs, and the report discards it.** SH-like support costs 141.9 s of
  239.96 s at 7 taxa x 2.07 Mbp (59%) and 8.5 s of 21.46 s at 4 taxa (40%) —
  but only 19.3 s of 545.6 s at 400 taxa x 100 kbp (3.5%), because it scales
  with alignment length where the ML search scales with taxon count. Few genomes
  over a megabase alignment is the pipeline's normal case, so the 59% is the
  figure that matters. `-nosupport` leaves topology and branch lengths
  byte-identical, and `draw_tree` labels leaves only, so the values never reach
  the page. Either drop them or render them — but `fasttree.newick` is a
  user-facing artifact too, so dropping them is a trade, not a free win.
- **Gap-filling per medium, parked 2026-09-03.** The idea that started
  `biosynthesis` was scoring genomes by growth on media standing in for
  ecological niches. That failed on measurement — none of the eleven Firmicute
  drafts grows on any defined medium, and a *B. subtilis* draft answers 29 of
  the 30 compounds and still returns zero, so one metabolite decides the score
  (measured 2026-09-10; a draft of *E. coli* does grow, so this is about the
  genome and not about CarveMe) — and the salvage of it was to gap-fill each
  model on each medium and
  score by *how many reactions had to be added*. Parked rather than dropped: it
  is a MILP per medium per genome on top of carving, so the cost is unmeasured
  and the SCIP presolver defect applies to it as well, and the CarveMe paper has
  already measured that the count confounds annotation quality with biology —
  gap-fill count against genome size, Pearson r = −0.29, P = 0.0055. Whoever
  picks it up should normalise against each genome's own minimum across media,
  and time one gap-fill on thylakoid before building anything.
- **Ecological niche media are not shipped, and that is not an oversight.** The
  compound sets for a gut, a rumen or a soil have no source that is both
  BiGG-mapped and checkable, and an invented one produces a figure nobody can
  audit. What ships instead is the panel and its compound classes; the
  aggregation that would name niches sits on top of per-compound calls, so it
  can be added when a citable source is.
- **Report sections** still to write from scratch: v2 never displayed fasttree or
  treecluster, so those are new rather than ports.
