# Installation

!!! warning "v3 is pre-release"
    The bioconda package is not published yet — the recipe is written and the
    code path works, but `conda install comparem2` still gives you v2. For now,
    install from git.

There are two ways to install CompareM2, and they differ in where the thirteen
analysis tools come from.

| | pixi, from git | conda *(not yet published)* |
| --- | --- | --- |
| the tools | one solved environment, all on PATH | deployed per rule on first run |
| first run | starts immediately | solves 14 environments |
| needs network at run time | only for databases | for databases and environments |
| good for | development, HPC, anything repeated | trying it once, or a shared install |

## Requirements

  - **Linux.** The analysis tools are `linux-64` only. On macOS you can run the
    unit tests and render reports, but not the pipeline.
  - **[pixi](https://pixi.prefix.dev/latest/#installation)** for the
    environment.
  - **Disk.** 1.58 GB of software, plus databases — see below. GTDB-Tk alone is
    60.8 GB.
  - **RAM.** GTDB-Tk's classify step is the peak; its own paper reports under
    55 GB for the v2 divide-and-conquer placement. Without GTDB-Tk, far less.

## Install

```bash
git clone https://github.com/cmkobel/CompareM2.git
cd CompareM2
pixi install
pixi run pytest        # 163 unit tests, no databases needed
```

!!! note "Moving the directory invalidates the environment"
    Conda bakes the absolute prefix into shebangs and RPATHs, so after moving
    the checkout you need `rm -rf .pixi && pixi install`.

## With conda, once published

The package will contain **the pipeline and none of the thirteen tools**.
Snakemake deploys each tool's own environment the first time it is needed:

```bash
conda install -c conda-forge -c bioconda comparem2
comparem2 *.fna --use-conda
```

`--use-conda` is the whole difference. Without it, CompareM2 expects the tools
on PATH and says so before doing any work:

```
not on PATH: seqkit (seqkit), mashtree (mashtree)
Either run inside an environment that has them (`pixi run cm2 ...`), or add
--use-conda to let Snakemake deploy each tool's own environment.
```

Environments go to `~/.comparem2/envs`, shared across runs, moved with
`--conda-prefix` or `$COMPAREM2_CONDA_PREFIX`. **Keep that path stable.**
Snakemake includes it in each environment's identity, so moving it re-solves all
fourteen — and re-fetches AMRFinder's database, which lives inside one of them.

Why the tools are not simply dependencies of the package: no single environment
can hold all thirteen. CheckM2 pins DIAMOND 2.1.x and current Bakta needs
2.2.x, so a package that listed them would have to leave one out.

## Environments

v3 solves **one** environment holding twelve of the thirteen tools. v2 shipped
25.

The exception is **CheckM2**, which pins DIAMOND 2.1.x while current Bakta needs
2.2.x — they cannot co-solve. CheckM2 gets its own environment, and because
Snakemake's shell does not inherit an interactive PATH, it needs an absolute
launcher:

```bash
pixi run cm2 *.fna \
  --isolated-launcher "$HOME/.pixi/bin/pixi run -e {tool}"
```

`{tool}` is substituted with the tool name. Without this, the CheckM2 step fails
with `command not found`.

## Databases

**Downloaded automatically.** Each database is a step in the workflow, so it is
fetched once, skipped if already present, and re-fetched if a previous attempt
was interrupted. Nothing needs to be placed by hand.

Sizes are measured from `content-length`, not estimated:

| Database | Size | Needed by |
|---|---|---|
| GTDB r232 | **60.8 GB** | `gtdbtk` |
| CheckM2 | 1.7 GB | `checkm2` |
| Bakta light | ~1.3 GB (documented, unmeasured) | `bakta` |
| AMRFinder | unmeasured | `amrfinder` |

They go under `databases/` unless you pass `-d`. Before running, CompareM2
prints what is actually missing — not a total that includes what you already
have:

```
2 assemblies, 5 tools
to download: amrfinder (1 of unknown size)
```

!!! warning "AMRFinder ignores `-d`"
    `amrfinder -u` refuses the `--database` option — *"only operates on the
    default database directory"* — so its data goes into the conda environment
    rather than your `--databases` directory. CompareM2 records that the update
    ran, but that record does not survive `pixi install` rebuilding the
    environment. Rebuild the environment and AMRFinder's database is fetched
    again.

**GTDB-Tk is 97% of the measured total**, and the release matters: 2.7 accepts
only r232 and refuses r226, which is also why the figure is 60.8 GB rather than
the 141.4 GB it was until r232 replaced FastANI's reference genomes with skani
sketches. It stays in the default path because it is the authoritative answer
to "what is this genome", but if you do not need taxonomy, leaving it out is
still the single biggest saving available:

```bash
pixi run cm2 *.fna --until seqkit checkm2 bakta amrfinder mlst \
                           mashtree treecluster skani panaroo snp-dists fasttree
```

Bakta uses the **light** database (~1.3 GB / 3.9 GB on disk) rather than v2's
full one (30 GB / 84 GB). That saves 29 GB for less specific functional
annotation, which a wide view can absorb — but note the Bakta paper's
annotation-quality figures are measured on the full database and do not
transfer.

## HPC

Snakemake's SLURM executor plugin is installed. Point CompareM2 at more cores
and let Snakemake submit:

```bash
pixi run cm2 *.fna --cores 64
```

Job submission is Snakemake's business, not CompareM2's — the generated
Snakefile lives at `<output>/.comparem2/Snakefile` and takes any Snakemake
profile you already use.

## Verification status

A clean `pixi install` says nothing about whether the pipeline runs. Both Bakta
and Panaroo once resolved to years-old builds that installed cleanly and crashed
on first use, which is why every tool now carries a minimum version.

`STATUS.md` in the repository tracks which command lines have actually been
**executed** on real genomes. All 13 now have. GTDB-Tk was the last and the
most instructive: its rule carried six defects that only running it could
show, including a database release the installed tool refuses.
