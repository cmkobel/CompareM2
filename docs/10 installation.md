
# Installation

## Install the bioconda package

<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/5b06b511-75c4-48cb-8ab8-f29b212ef6df">

It is highly recommended that you have [Apptainer](https://Apptainer.org/docs/user/main/quick_start.html#installation-request) on your system as it makes Assemblycomparator2 able to use a compressed Docker-image that speeds up installation significantly.

<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/c9d15678-b4a7-42be-b0de-b649479f6d74">

First, you need to install a Conda or Mamba package manager.
The recommended choice is [Miniforge](https://github.com/conda-forge/miniforge#install) which not only provides the required Python and Conda commands, 
but also includes Mamba - an extremely fast and robust replacement for the Conda package manager which is highly recommended.

!!! info
    In case you don't use Miniforge you can always install Mamba into any other Conda-based Python distribution with:

    ```
    conda install -n base -c conda-forge mamba
    ```


<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/6bc39697-7e90-49a0-a44e-64820f2c1024">

Finally, Assemblycomparator2 can be installed into its own environment with the correct channels like so:

```bash

mamba create -c conda-forge -c bioconda -n asscom2 assemblycomparator2

```

Installing into isolated environments is best practice in order to avoid side effects with other packages.


!!! note
    If you want to develop new rules in the Assemblycomparator2 pipeline, you should consider following [the development version installation instructions](https://github.com/cmkobel/assemblycomparator2/blob/master/readme-development.md). The development version contains the full git repository and is purely conda-based so you can affect the next version of the Apptainer-compatible Docker image. 


## Optionally: Testing the installation

Now you will be able to run asscom2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run asscom2 it will show the message "Pulling singularity image docker://cmkobel/assemblycomparator2." This might take some time depending on your network bandwidth as it downloads a +4GB Docker image that contains all the conda environments needed for each analysis.

```bash

# Activate the newly created conda environment containing the asscom2 launcher.
mamba activate asscom2

# First, create an empty directory and enter.
mkdir test_ac2_install
cd test_ac2_install

# Copy some test metagenomic assemblies from the test directory.
cp $CONDA_PREFIX/assemblycomparator2/tests/MAGs/*.fasta .

# Should take about a minute to complete the "fast" pseudo-rule.
asscom2 --until fast

# You can then investigate the report document that has been generated.
# open results_ac2/report_test_ac2_install.html

# Downloads all databases (~ 200 GB).
asscom2 --until downloads

# Run the full pipeline (~ 1 cpu-hour per genome).
asscom2
```



## Advanced configuration

### Shared database

If you are working on a shared computational resource like a laboratory workstation or a HPC you might want to share a database directory so that each user will not have to redundantly download each database. To set this up, the first user must decide on a directory and set reading and writing permissions for the group of users that should be able to use the database. Writing permissions are necessary for the "database representative" flags that snakemake uses to keep track of the presence of the databases. Setting this custom path is a matter of defining the "ASSCOM2_DATABASES" environment variable. You can put this into your ~/.bashrc or execute the command before using asscom2.

```bash
export ASSCOM2_DATABASES="/absolute/path/to/shared_databases/asscom2_v2.5.8+"
```

### HPC profiles for Snakemake

 If you have experience with snakemake and are working on a high performance computing cluster (HPC), you can modify and use the cluster configuration profiles in the "profiles/" directory. You can define the use of one of these profiles by setting the "ASSCOM2_PROFILE" environment variable. You can put this into your ~/.bashrc or execute the command before using asscom2. You can read more about snakemake profiles [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) or browse more default profiles [here](https://github.com/snakemake-profiles).

```bash
export ASSCOM2_PROFILE=${ASSCOM2_BASE}/profiles/apptainer/slurm-sigma2-saga
```






