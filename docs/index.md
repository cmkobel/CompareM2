# CompareM2

CompareM2 takes microbial genome assemblies — isolates or MAGs, from any
sequencing technology — and produces one portable HTML report comparing them.

!!! note "Looking for something else?"
    The paper describes v2, documented at
    [comparem2.readthedocs.io/en/stable](https://comparem2.readthedocs.io/en/stable/).

    For the original CompareM (AAI and codon usage), see
    [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM).

## Built for many genomes at once

Fourteen analyses run across a whole set of assemblies and land in one report,
and the cost scales near-linearly with how many you give it — which is the
benchmark result the design is built on. Comparing a hundred genomes is the
case it is good at, and for most questions the per-genome depth is what you
need. Where you want a deep functional or metabolic dive into a single genome,
[DRAM2](https://github.com/WrightonLabCSU/DRAM) and
[nf-core/funcscan](https://nf-co.re/funcscan) complement it well.

Narrowing the scope is what pays for everything else: **14 tools instead of
30+, two conda environments instead of 25, 7.7 GB of software.**

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
  - [Installation](10 installation.md) — conda, pixi, databases, HPC
  - [Usage](20 usage.md) — the CLI, the TUI, passthrough parameters
  - [What analyses does it do](30 what analyses does it do.md) — the 14 tools
  - [Citing](99 citation.md) — CompareM2 and every tool it runs
