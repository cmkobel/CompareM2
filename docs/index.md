# CompareM2

CompareM2 takes microbial genome assemblies — isolates or MAGs, from any
sequencing technology — and produces one portable HTML report comparing them.

!!! note
    Looking for the original CompareM (AAI and codon usage)? See
    [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM).

## Built for many genomes at once

Fourteen analyses run across the whole set and land in one report, and the cost
scales near-linearly with how many assemblies you give it — the benchmark
result the design is built on. A hundred genomes is an ordinary input, not a
stress test.

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
  - **Metabolism** — draft genome-scale metabolic models (CarveMe), and what
    each one can build for itself: 32 amino acids, vitamins and cofactors, de
    novo or acquired

See [what analyses does it do](30 what analyses does it do.md) for the full
reference, generated from the tool specs themselves.

## An interpretative report

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

Installed it already and want to see a report before committing any genomes of
your own? `comparem2 --demo` runs four analyses over six bundled *Enterococcus
faecium* plasmids and downloads nothing.

  - [Quick start](05 quick start.md) — install and run
  - [Installation](10 installation.md) — conda, pixi, databases, HPC
  - [Usage](20 usage.md) — the CLI, the TUI, passthrough parameters
  - [What analyses does it do](30 what analyses does it do.md) — the 14 tools
  - [Citing](99 citation.md) — CompareM2 and every tool it runs
