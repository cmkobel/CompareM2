# CompareM2 v3 — design decisions

Working notes for the v3 rewrite. Decisions are dated; if one is reversed, edit
the entry and say why rather than deleting it.

## What survives from v2

Only the philosophy: **easy to install, easy to run, easy to interpret** — a
straightforward view into a set of assemblies. Everything else is a blank slate.

## Decisions

### 2026-09-01 — Lives on branch `v3` in `cmkobel/comparem2`
Becomes CompareM2 v3 rather than a separate project, so the name, citation and
docs domain carry over. v2 code is left in place on this branch for reference
and will be removed in a deliberate, separate step once v3 stands on its own.

### 2026-09-01 — Python only; R is dropped
CLI, TUI and report all share one runtime. Removes R, pandoc, rmarkdown,
tidyverse and r-clusterProfiler from the install path — the heaviest
non-database dependency in v2. The 16 existing `.rmd` report sections are not
ported; they are rewritten.

### 2026-09-01 — Snakemake is kept, as an executor rather than a framework
Reversed within the day. The first call was to drop Snakemake on the grounds
that a single-digit tool count makes a DAG engine overkill. That reasoning was
wrong: the DAG was never the hard part. The expected deployment is a workstation
or an HPC head node with jobs submitted to SLURM when available, and *cluster
execution* is the hard part — sbatch generation, `afterok` dependency chains,
`squeue`/`sacct` polling, partial-failure recovery, retries. Snakemake's SLURM
executor plugin already solves it. The install-weight argument for dropping it
also fails on the numbers: Snakemake is tens of MB against 141 GB of databases.

Snakemake is a Python package, so this does not conflict with dropping R.

Shape: the `Tool` specs in `src/comparem2/tools.py` are the single source of
truth, and Snakemake rules are **generated** from them. Adding a tool stays
"add one spec"; Snakemake owns DAG resolution, resumability, retries and
SLURM submission. CLI, TUI and report all read the same specs.

Cost, accepted knowingly: **the TUI gets harder.** Snakemake owns the event loop
and its logging, so a live keyboard-driven interface has to drive it through
`snakemake.api` with a custom log handler, and depends on that log event schema.

### 2026-09-01 — Install weight is a first-class constraint
v2 shipped 25 conda environments and 7 databases. GTDB-Tk's r226 full package
alone is 141.4 GB (measured from `content-length`, 2026-09-01) against CheckM2's
1.7 GB — a factor of 82. Any tool that would be added to v3 declares its
database size in bytes, and the total is shown to the user before download.

### 2026-09-01 — Wide, not deep
v3 is a quick wide view over many assemblies, not a deep view into each one.
Deep functional and metabolic analysis is where DRAM2 and nf-core/funcscan are
already strong; triage across a large set is unoccupied, and it is the direction
the v2 benchmark result (scaling with input count) actually supports.

Consequence: gtdbtk, antismash, gapseq, interproscan and eggnog are out of the
default path. That is what makes the next decision possible.

### 2026-09-01 — One environment, plus CheckM2 on its own (revised same day)
**The first version of this decision was wrong, and the way it was wrong is
worth keeping.** The original solve below succeeded, so "one environment" was
recorded as settled. It only succeeded because the solver quietly selected
**bakta 1.8.1** (2023) — the newest Bakta that could co-exist with CheckM2.
Bakta 1.8.1 then crashes at runtime against pyrodigal 3.x, which renamed
`OrfFinder` to `GeneFinder`. **A solve is not a working environment**, and only
running the tool revealed it.

Isolated afterwards, on `linux-64`:

| Combination | Result |
| ----------- | ------ |
| bakta ≥1.10 alone | solves — 1.12.1, diamond 2.2.5 |
| bakta ≥1.10 + checkm2 | **conflict** |
| bakta ≥1.10 + the other thirteen tools | solves |

CheckM2 pins DIAMOND 2.1.x; current Bakta needs 2.2.x. CheckM2 is too valuable
to drop (1.7 GB for completeness and contamination), so it is the one tool
marked `isolated=True` and given its own environment. Thirteen tools share the
main environment.

Rule: `isolated` is an exception that must carry its reason in the spec.
v2 reached 25 environments by treating it as the default.

### 2026-09-01 — The original (superseded) single-environment measurement
Measured with pixi against `linux-64`, not assumed:

| Environment                | Packages | Download |
| -------------------------- | -------: | -------: |
| core 8 tools               |      576 |  0.76 GB |
| core + bakta + snakemake   |      830 |  0.83 GB |
| core + prokka + snakemake  |      840 |  1.60 GB |

(core = seqkit, checkm2, mashtree, treecluster, mlst, ncbi-amrfinderplus,
sylph, skani.) All three solve with no conflicts, against v2's 25 environments.
This works *because* the wide decision removed the historical conflict sources;
re-test before adding any of them back. Note prokka costs roughly twice bakta in
package weight — it pulls perl and BLAST.

### 2026-09-01 — The tool set: 14 in, 8 out
Selected tool by tool. v2 had 20; v3 has 14, of which 2 are new.

**Kept** — seqkit, checkm2, gtdbtk, sylph, bakta, amrfinder, mlst, mashtree,
treecluster, skani, panaroo, snp-dists, fasttree, carveme.

**Dropped** — assembly-stats (seqkit covers it), prokka (bakta only; prokka cost
0.84 GB more in packages, measured), eggnog, dbcan, interproscan, gapseq,
iqtree, clusterProfiler (an R package, and R is gone).

**New** — sylph (fast provisional taxonomy while GTDB-Tk runs) and skani
(all-against-all ANI).

Two decisions worth their reasoning:
- **carveme in, gapseq out.** Genome-scale metabolic models are the one
  capability Bactopia, funcscan and DRAM2 all lack. CarveMe keeps it at minutes
  rather than gapseq's hours. It solves with open-source SCIP — no CPLEX
  licence — but it was commented out of v2 (`Snakefile:332`), so it needs
  verifying that it still runs.
- **antismash dropped.** It was selected, then dropped when it turned out to
  break the single environment: it pins `biopython 1.78` and `diamond 2.1.11`
  against the newer numpy/biopython that checkm2, bakta and gtdbtk require.
  One environment was judged worth more than one BGC caller, especially as
  funcscan ships four.

Final environment, measured for `linux-64`: **1142 packages, 1.58 GB, no
conflicts.** Against v2's 25 environments.

### 2026-09-01 — Bakta light, not full
v2 downloads `--type full` (30 GB compressed, 84 GB on disk). v3 uses `light`
(1.3 GB / 3.9 GB, Bakta's documented figures). Saves 29 GB for less specific
annotation, which a wide view can absorb.

### 2026-09-01 — Vertical slice runs end to end
seqkit → mashtree → treecluster on v2's four *E. faecium* test genomes,
on thylakoid (Linux, 24 cores). 7/7 steps, report rendered. This closes the
central architectural question: **Snakemake rules can be generated from
declarative specs**, including wildcards, dependencies and stdout redirection.

Sanity checks that the result is real, not merely present: mashtree puts
`116_2` and `116_2_duplicate` at distance 0.00000, and seqkit derives identical
statistics for them (2,712,340 bp, 4 contigs, 38.0% GC).

Three bugs the slice caught that no amount of reading would have:
1. **Bare `{sample}` in a shell block is a runtime NameError.** Snakemake
   resolves wildcards in `input:`/`output:` but exposes them as
   `wildcards.sample` inside `shell:`. The Snakefile parsed and built the right
   DAG, then failed on execution.
2. **`TreeCluster.py` requires `-t/--threshold`.** The drafted command omitted
   it. v2's defaults (`max_clade`, `0.05`) were adopted.
3. **Sample names can contain spaces.** v2's own test set ships
   `116_2 duplicate.fna`, which produces a silently broken wildcard.

### 2026-09-01 — Passthrough parameters survive from v2
v2's `set_<tool>--<flag>: <value>` is now `--set tool--flag=value`, backed by a
`params` field on each `Tool`. Defaults are carried over from v2's
`config.yaml` verbatim, so v3 reproduces v2's behaviour unless told otherwise.

### 2026-09-01 — Unit tests, which v2 never had
v2 validated only by running the pipeline in CI. That is the wrong instrument
for a generator: a wrong wildcard yields a Snakefile that parses cleanly and
builds the wrong DAG. 26 tests cover name sanitisation, dependency closure,
rule generation and report rendering. `pixi run pytest tests/unit -q`.

## Command-line verification status

Drafted commands are worthless until executed. Verified means: ran on v2's four
*E. faecium* test genomes on thylakoid and produced correct-looking output.

| Tool        | Status     | Note                                             |
| ----------- | ---------- | ------------------------------------------------ |
| seqkit      | verified   | per-contig lengths and GC                         |
| mashtree    | verified   | duplicates at distance 0.00000, as they should be |
| treecluster | verified   | needed `--threshold`; v2 defaults adopted         |
| skani       | verified   | duplicates at 100.00% ANI                         |
| mlst        | verified   | ST32 / ST117 / ST78, correct for *E. faecium*     |
| checkm2     | verified   | 100% complete; `--database_path` wants the .dmnd **file** |
| bakta       | verified   | needed `--force` and a pin to >=1.10; db must be v6 |
| carveme     | verified   | **it does run** — ~1200 reactions/genome, but ~9 min each |
| amrfinder   | verified   | protein-only, as v2; 7/11/8 genes, incl. 5 glycopeptide in E8202 |
| panaroo     | verified   | needed pin >=1.5; 3780 clusters, 2091 core        |
| snp-dists   | verified   | 0 SNPs between the duplicate pair                 |
| fasttree    | verified   | duplicates at 0.0 branch length                   |
| gtdbtk      | unverified | needs 141.4 GB; thylakoid has 180 GB free         |
| sylph       | unverified | GTDB sketch release still unchosen                |

**12 of 14 verified.** The two outstanding both concern taxonomy: GTDB-Tk needs
the 141.4 GB download, and sylph needs a database release to be chosen.

Cross-checks used throughout: v2's test set contains `116_2.fna` and
`116_2 duplicate.fna`, the same genome twice. Any tool that treats them
differently is wrong. Every verified tool agrees — 0.00000 mash distance,
100.00% ANI, 0 SNPs, identical CDS counts, identical model sizes.

## Where the install weight actually sits

| Database   | Download   | Note                             |
| ---------- | ---------: | -------------------------------- |
| GTDB-Tk    | 141.4 GB   | measured; 91% of the total       |
| CheckM2    |    1.7 GB  | measured                         |
| Bakta light|    1.3 GB  | 3.2 GB on disk, measured          |
| AMRFinder  |  unmeasured| may be free — Bakta's light DB ships an `amrfinderplus-db` |
| sylph GTDB |  unmeasured| release not yet chosen           |
| *software* |    1.58 GB | measured, one environment        |

The environment is slim. The data is not, and GTDB-Tk is the whole reason.
Making it opt-in was offered twice and declined twice; it stays in the default
path. If "easy to install" ever has to be defended, this is the line item.

### 2026-09-02 — Pin a minimum version for every tool
Two tools were silently resolved to years-old builds that install cleanly and
crash on first use:

| Tool | Unpinned resolved to | Breaks on |
| ---- | -------------------- | --------- |
| bakta | 1.8.1 (2023) | `pyrodigal.OrfFinder`, renamed `GeneFinder` in pyrodigal 3.x |
| panaroo | 1.1.2 (2020) | `Bio.Alphabet`, removed in Biopython 1.78 |

Neither is a solver failure — the solver did what it was asked. Bioconda keeps
old builds forever, and an unconstrained `tool = "*"` lets the solver satisfy
some *other* package's constraint by reaching back years. **`pixi install`
succeeding says nothing about whether the pipeline works.**

Hence: pin a minimum version in each spec's `conda` field, and treat the
verification table above as tracking *execution*, never installation.

### 2026-09-02 — CarveMe is real, and it is the slow one
CarveMe was commented out in v2 and its status was an open question. It runs:
~4 MB SBML model per genome, all four test genomes. But it takes roughly nine
minutes per genome single-threaded, which makes it the slowest per-genome tool
in the set now that InterProScan is gone. For a wide view over hundreds of
assemblies that matters, and it should probably be opt-in rather than default.

### 2026-09-02 — Working checkout lives on thylakoid
`/evo/postdoc/cm2v3` (reachable as `~/postdoc/cm2v3`). 24 cores, 125 GB RAM.
Databases are in `databases/`; `pixi run cm2 …` is the entry point, and
isolated tools need `--isolated-launcher "/home/thylakoid/.pixi/bin/pixi run
-e {tool}"` with an absolute path, because Snakemake's shell does not inherit
an interactive PATH.

Moving a pixi project invalidates its environments — conda bakes the absolute
prefix into shebangs and RPATHs — so `rm -rf .pixi && pixi install` is required
after any move. Results directories also hold absolute symlinks to the inputs;
`canonicalise()` now repairs dangling ones rather than raising.

GenomeDK was the alternative and is not reachable non-interactively from here
(`Permission denied (publickey,keyboard-interactive)`).

## Open questions
- **Taxonomy.** The single largest install cost. GTDB-Tk is authoritative and
  141.4 GB; sylph/skani/mash-against-GTDB are fast and small but approximate.
  Candidates not yet measured.
- **Every command line in `catalogue.py` is unverified.** They were written
  against documented interfaces, never executed — the tools are linux-64 and
  the selection was made on macOS. This is the largest single risk in the
  rewrite and should be the next thing closed.
- **Three database sizes are unmeasured**: bakta-light (documented as 1.3 GB but
  not measured), amrfinder, and sylph's GTDB sketch. Measure before any total is
  shown to a user.
- **sylph's GTDB release** is unchosen, and it should match GTDB-Tk's r226 or
  the two tools will disagree about taxonomy in confusing ways.
- **CarveMe needs verifying**: that it runs at all (it was disabled in v2), and
  whether DIAMOND must be supplied explicitly — it is absent from CarveMe's
  dependency closure, and v2 shipped a separate `diamond` environment.
- **GTDB-Tk's summary output.** It writes `bac120` and `ar53` summaries
  separately; the declared `gtdbtk.summary.tsv` requires a concatenation step
  that does not exist yet.
- ~~**TUI against Snakemake.**~~ **Resolved 2026-09-01 — no fork needed.**
  Snakemake 9.26.1 ships a logger plugin system (`LogHandlerBase`, `LogEvent`,
  `LoggerPluginRegistry`), and events can also be captured in-process by
  attaching a `logging.Handler` to the `snakemake` logger while driving
  `SnakemakeApi`. A spike captured 19 events on a two-job run, with structured
  fields rather than scraped text:

  | Event | Carries |
  | ----- | ------- |
  | `run_info` | per-rule job counts and total |
  | `job_info` | `jobid`, `rule_name` |
  | `job_started` / `job_finished` | job lifecycle |
  | `progress` | `done`, `total` |
  | `shellcmd` | the exact command |
  | `job_error` / `error` | failures |

  That is everything a live keyboard-driven interface needs. The offer to ship
  a modified Snakemake stays unused, which is the better outcome — a fork would
  have to be rebased forever.
- **Report sections** for all 14 tools have to be written from scratch in
  Python. Note v2 never displayed antismash, interproscan, iqtree, fasttree or
  treecluster — so fasttree and treecluster need genuinely new sections, not
  ports.
