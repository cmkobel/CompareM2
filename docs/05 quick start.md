# Quick start

Linux only. The analysis tools are `linux-64`.

With [pixi](https://pixi.prefix.dev/latest/#installation):

```bash
pixi global install --channel conda-forge --channel bioconda comparem2
```

or with conda:

```bash
conda install -c conda-forge -c bioconda comparem2
```

Either puts `comparem2` on your `PATH`, along with `cm2` as a shorter alias for
it. For workspace-scoped installs, see [installation](10 installation.md).

!!! info "Working on CompareM2 rather than using it?"
    Install from git with pixi instead (see
    [installation](10 installation.md)) and read every `comparem2` below as
    `pixi run comparem2`.

## 1) Run it

On your own assemblies:

```bash
comparem2 *.fna
```

Before fetching anything, CompareM2 says what it is about to do and what it
will cost:

```
4 assemblies, 14 tools
to download: checkm2, gtdb, bakta-light, amrfinder (62.5 GB + 2 of unknown size) -> /home/you/.comparem2/databases
tool environments: 6 in /home/you/.comparem2/envs (none built yet)
```

Of that 62.5 GB, GTDB-Tk accounts for 60.8 GB. If you do not need taxonomic
assignment, name the analyses you do want and the download drops to roughly
3 GB; see [running a subset](20 usage.md#running-a-subset). The figure printed
there reads 1.7 GB, because Bakta's and AMRFinder's databases publish no size.
They cost 1.3 GB and 241 MB.

The first run is also quiet for about a minute while Snakemake builds those six
tool environments. That is expected, and `comparem2 --setup` does it in advance.

!!! tip "No genomes to hand?"
    `comparem2 --demo` needs none, and no databases either. Six *Enterococcus
    faecium* plasmids ship inside the package, and it runs the four analyses
    that need no database. See
    [the bundled demo](20 usage.md#the-bundled-demo).

## 2) Read the report

`results_comparem2/report.html` is one self-contained file with no external
assets, so it survives being emailed or copied off a cluster.

Each section carries a **"What this is, and how to read it"** block, collapsed
by default. Open it before drawing a conclusion from the numbers. It says what
the columns mean, what the tool's own error is, and what the result cannot
show.

[An example report](07 an example report.md) is a real run of all fourteen
tools over eight genomes, with a walk through what to look at first.

---

[Installation](10 installation.md) covers databases, the tool environments and
HPC. [Usage](20 usage.md) covers the full CLI, the TUI and passthrough
parameters.
