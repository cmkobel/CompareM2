
# Installation

## Install with pixi (recommended)

<img width="150" align="right" src="https://github.com/cmkobel/comparem2/assets/5913696/5b06b511-75c4-48cb-8ab8-f29b212ef6df">

It is recommended to have [Apptainer](https://Apptainer.org/docs/user/main/quick_start.html#installation-request) on your system. It allows CompareM2 to use a pre-built Docker image, which speeds up installation significantly.

<img width="150" align="right" src="https://raw.githubusercontent.com/cmkobel/CompareM2/refs/heads/master/resources/pixi-logo.svg">

First, install [pixi](https://pixi.prefix.dev/latest/#installation) — a fast package manager for conda-forge and bioconda packages.

Then install CompareM2 globally:

```bash
pixi global install -c conda-forge -c bioconda comparem2
```

This creates an isolated environment and makes the `comparem2` command available on your PATH.

!!! note
    If you want to develop new rules in the CompareM2 pipeline, you should consider following [the development version installation instructions](https://github.com/cmkobel/comparem2/blob/master/readme-development.md). The development version contains the full git repository and is purely conda-based so you can affect the next version of the Apptainer-compatible Docker image.

## Alternative: Install with mamba

If you prefer conda/mamba, install CompareM2 with [Miniforge](https://github.com/conda-forge/miniforge#install):

```bash
mamba create -c conda-forge -c bioconda -n comparem2 comparem2
```

Activate the environment before each use:

```bash
mamba activate comparem2
```


## Testing the installation

After installation, verify everything works using the bundled test data:

```bash
# Create a fresh working directory
mkdir test_comparem2_install
cd test_comparem2_install

# Extract test genomes from the installation
unzip $CONDA_PREFIX/share/comparem2-latest/tests/E._faecium/fna.zip

# Run the fast pseudo-rule (should complete in about a minute)
comparem2 --until fast

# Open the generated report
# open results_comparem2/report_test_comparem2_install.html

# Download all databases (~200 GB)
comparem2 --until downloads

# Run the full pipeline (~1 cpu-hour per genome)
comparem2
```


## Advanced configuration

### Shared database directory

On shared systems (lab workstations, HPC clusters), multiple users can share a single database directory to avoid redundant downloads. Set the `COMPAREM2_DATABASES` environment variable to a shared path with group read/write permissions:

```bash
export COMPAREM2_DATABASES="/absolute/path/to/shared_databases/comparem2_db"
```

Add this to your `~/.bashrc` to make it persistent. System administrators can set it globally in `/etc/bash.bashrc`.

### HPC profiles for Snakemake

If you work on a high-performance computing (HPC) cluster, you can use or customize the cluster profiles in the `profiles/` directory. Set the `COMPAREM2_PROFILE` environment variable to point to your profile:

```bash
export COMPAREM2_PROFILE=${COMPAREM2_BASE}/profiles/apptainer/slurm-sigma2-saga
```

See the [Snakemake profiles documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) or browse [community profiles](https://github.com/snakemake-profiles).


{!resources/footer.md!}
