# CompareM2

CompareM2 takes microbial genome assemblies — isolates or MAGs, from any
sequencing technology — and produces one portable HTML report comparing them.

!!! warning "These docs describe v3, which is pre-release"
    v3 is a rewrite and there is no released package for it yet — install from
    git, see [installation](10 installation.md). The version the paper
    describes is v2, documented at
    [comparem2.readthedocs.io/en/stable](https://comparem2.readthedocs.io/en/stable/).

!!! note
    Looking for the original CompareM (AAI and codon usage)? See
    [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM).

## A wide view, not a deep one

v3 is triage across many assemblies rather than a deep dive into each. Deep
functional and metabolic analysis is where DRAM2 and nf-core/funcscan are
already strong; a fast, wide sweep over a large set is not well served, and it
is what CompareM2's benchmark result — near-linear scaling with input count —
actually supports.

Narrowing the scope is what pays for everything else: **14 tools instead of
30+, one conda environment instead of 25, 1.58 GB of software.**

## What it does

  - **Quality** — contig lengths, GC and N50 (SeqKit); completeness and
    contamination (CheckM2)
  - **Taxonomy** — GTDB-Tk against the Genome Taxonomy Database
  - **Annotation** — Bakta
  - **Screening** — antimicrobial resistance and virulence (AMRFinderPlus),
    sequence typing (MLST)
  - **Relatedness** — alignment-free tree (Mashtree), clusters (TreeCluster),
    all-against-all ANI (skani)
  - **Pangenome** — core and accessory gene content (Panaroo), SNP distances,
    core-genome tree (FastTree)
  - **Metabolism** — draft genome-scale metabolic models (CarveMe)

See [what analyses does it do](30 what analyses does it do.md) for the full
reference, generated from the tool specs themselves.

## The report explains itself

Every section says what the tool does, how to read the specific columns on
screen, and what the result *cannot* tell you — with the numbers quoted from
each tool's own paper and checked against it. The report ends with a *Methods
and citations* list covering exactly the tools that ran, which is the thing you
paste into a manuscript.

That matters more than it sounds. A completeness of 92% and one of 94% are not
meaningfully different — CheckM2's own mean absolute error is 2.1±2.9% — and a
report that prints both without saying so invites a conclusion the data will
not support.

## Where to start

  - [Quick start](05 quick start.md) — install and run
  - [Installation](10 installation.md) — pixi, databases, HPC
  - [Usage](20 usage.md) — the CLI, the TUI, passthrough parameters
  - [What analyses does it do](30 what analyses does it do.md) — the 14 tools
  - [Citing](99 citation.md) — CompareM2 and every tool it runs
