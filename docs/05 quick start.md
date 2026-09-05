# Quick start

Linux only — the analysis tools are `linux-64`.

With [pixi](https://pixi.prefix.dev/latest/#installation):

```bash
pixi global install --channel conda-forge --channel bioconda comparem2
```

or with conda:

```bash
conda install -c conda-forge -c bioconda comparem2
```

Either puts `comparem2` and `cm2` on your `PATH`. See
[installation](10 installation.md) for workspace-scoped installs.

!!! info "Working on CompareM2 rather than using it?"
    Install from git with pixi instead — see
    [installation](10 installation.md) — and read every `cm2` below as
    `pixi run cm2`.

## 1) Run it

On your own assemblies:

```bash
cm2 *.fna
```

Before fetching anything, CompareM2 says what it is about to do and what it
will cost:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> /home/you/.comparem2/databases
tool environments: 2 in /home/you/.comparem2/envs (none built yet)
```

**60.8 GB of that 62.5 GB is GTDB-Tk alone.** If you do not need taxonomic
assignment, name the analyses you do want instead and the download falls under
2 GB — see [running a subset](20 usage.md#running-a-subset).

The first run is then quiet for about a minute while Snakemake builds those two
tool environments. That is expected, and `cm2 --setup` does it in advance.

!!! tip "No genomes to hand?"
    `cm2 --demo` needs none, and no databases either. Six *Enterococcus
    faecium* plasmids ship inside the package, and one of them is passed a
    second time under another name — so seven inputs, six distinct genomes. It
    runs the four analyses that need no database, and the duplicated pair is
    what you check the report against: it must come out at 0.00000 mash
    distance and 100.00% ANI.

## 2) Read the report

`results_comparem2/report.html` — one self-contained file with no external
assets, so it survives being emailed or copied off a cluster.

Each section carries a **"What this is, and how to read it"** block, collapsed
by default. Open it before drawing a conclusion from the numbers: it says what
the columns mean, what the tool's own error is, and what the result cannot show.

---

[Installation](10 installation.md) covers databases, the tool environments and
HPC. [Usage](20 usage.md) covers the full CLI, the TUI and passthrough
parameters.
