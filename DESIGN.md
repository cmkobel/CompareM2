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

## Open questions

- **Wide vs deep.** "Quick wide view over many assemblies" and "deep view into
  each assembly" imply different tool sets. Unresolved; governs tool selection.
- **Taxonomy.** The single largest install cost. GTDB-Tk is authoritative and
  141.4 GB; sylph/skani/mash-against-GTDB are fast and small but approximate.
  Candidates not yet measured.
- **Tool set.** Selected manually by Carl. Not yet chosen.
- **One conda environment or several?** v2 solved 25 environments on first run,
  which dominated install time. A slim set may fit in one, killing that cost —
  but only if the chosen tools' dependencies co-solve. Verify once the tool set
  exists; a single container is the fallback.
- **TUI against Snakemake.** Needs a spike: drive `snakemake.api` with a custom
  log handler and confirm the event stream is rich enough for live progress.
