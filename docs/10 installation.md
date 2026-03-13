
# Installation

## Install with pixi (recommended)

<img width="150" align="right" src="https://github.com/cmkobel/comparem2/assets/5913696/5b06b511-75c4-48cb-8ab8-f29b212ef6df">

It is recommended that you have [Apptainer](https://Apptainer.org/docs/user/main/quick_start.html#installation-request) on your system as it makes CompareM2 able to use a compressed Docker-image that speeds up installation significantly.

<img width="150" align="right" src="https://github.com/cmkobel/comparem2/assets/5913696/c9d15678-b4a7-42be-b0de-b649479f6d74">

First, you need to install [pixi](https://pixi.sh) — a fast package manager for conda-forge and bioconda packages.

Then, CompareM2 can be installed globally in a single step:

```bash

pixi global install -c conda-forge -c bioconda comparem2

```

This creates an isolated environment and makes the `comparem2` command available on your PATH.

!!! note
    If you want to develop new rules in the CompareM2 pipeline, you should consider following [the development version installation instructions](https://github.com/cmkobel/comparem2/blob/master/readme-development.md). The development version contains the full git repository and is purely conda-based so you can affect the next version of the Apptainer-compatible Docker image.

## Alternative: Install with mamba

If you prefer using conda/mamba, CompareM2 can also be installed with [Miniforge](https://github.com/conda-forge/miniforge#install):

```bash

mamba create -c conda-forge -c bioconda -n comparem2 comparem2

```

You will then need to activate the environment before using CompareM2:

```bash

mamba activate comparem2

```


## Optionally: Testing the installation

Now you will be able to run CompareM2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run CompareM2 it will show the message "Pulling singularity image docker://cmkobel/comparem2." This might take some time depending on your network bandwidth as it downloads a +4GB Docker image that contains all the conda environments needed for each analysis.

```bash

# First, create an empty directory and enter.
mkdir test_comparem2_install
cd test_comparem2_install

# Copy some test metagenomic assemblies from the test directory.
# cp $CONDA_PREFIX/share/comparem2-*/tests/E._faecium/*.fna . # Until v2.7.1 you must copy the .fna files
# unzip $CONDA_PREFIX/share/comparem2-*/tests/E._faecium/fna.zip # From 2.8.1
unzip $CONDA_PREFIX/share/comparem2-latest/tests/E._faecium/fna.zip # From 2.12.1

# Should take about a minute to complete the "fast" pseudo-rule.
comparem2 --until fast

# You can then investigate the report document that has been generated.
# open results_comparem2/report_test_comparem2_install.html

# Downloads all databases (~ 200 GB).
comparem2 --until downloads

# Run the full pipeline (~ 1 cpu-hour per genome).
comparem2
```



## Advanced configuration

### Shared database

If you are working on a shared computational resource like a laboratory workstation or a HPC you might want to share a database directory so that each user will not have to redundantly download each database. To set this up, the first user must decide on a directory and set reading and writing permissions for the group of users that should be able to use the database. Writing permissions are necessary for the "database representative" flags that snakemake uses to keep track of the presence of the databases. Setting this custom path is a matter of defining the "COMPAREM2_DATABASES" environment variable. You can put this into your ~/.bashrc or execute the command before using CompareM2.

```bash
export COMPAREM2_DATABASES="/absolute/path/to/shared_databases/comparem2_db"
```

As a system administrator you may wish to set the "COMPAREM2_DATABASES" variable globally for all users. On common debian/ubuntu systems using bash, this can be done by putting the same export statement into /etc/bash.bashrc

### HPC profiles for Snakemake

 If you have experience with snakemake and are working on a high performance computing cluster (HPC), you can modify and use the cluster configuration profiles in the "profiles/" directory. You can define the use of one of these profiles by setting the "COMPAREM2_PROFILE" environment variable. You can put this into your ~/.bashrc or execute the command before using CompareM2. You can read more about snakemake profiles [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) or browse more default profiles [here](https://github.com/snakemake-profiles).

```bash
export COMPAREM2_PROFILE=${COMPAREM2_BASE}/profiles/apptainer/slurm-sigma2-saga
```






{!resources/footer.md!}
