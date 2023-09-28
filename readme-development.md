# Development information

This file describes how to install assemblycomparator2 for development purposes. Since you will be defining and configuring a conda environment to add to the apptainer-compatible dockerfile, you must run everything through conda. Once you make a pull request and it is accepted, your conda environment will be added to the official docker image which is then published with a new version number.

The installation process is more or less the same as with the normal installation, but you must make sure to use the conda configuration. 

Below is a copy of the normal installation procedure, but where the changes that have been made for the development version (purely conda based) is annotated.

## Installation of assemblycomparator2 on Linux

assemblycomparator2 can be installed by downloading the code and setting up an alias in your user profile that let's you launch the pipeline from any directory on your machine.

The only requisites for running assemblycomparator2 is:
  - [*miniconda*](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) package manager
  - *git* distributed version control (can be installed with conda by typing `conda install -c anaconda git`)
  - ~~*apptainer* container-virtualizer~~


#### 0) Prerequisites

First, check that you have the prerequisites available on your system. Open a terminal and paste these commands. 

```bash
which conda && conda --version
which git && git --version
# which apptainer && apptainer --version # Not necessary for the conda-based development version.
```

#### 1) Download pipeline and set up the launcher environment

Then download the assemblycomparator2 pipeline and set up an alias in your profile (.bashrc on most linux distributions). Proposed installation directory is in your home directory (\~).

```bash
cd ~ # Enter the directory where you want to install assemblycomparator2.
git clone https://github.com/cmkobel/assemblycomparator2.git asscom2
cd asscom2
conda env create --name asscom2_launcher --file environment.yaml # Installs snakemake and mamba in an environment named "asscom2_launcher".

```


#### 2) Alias

Finally, define the alias that will be used to launch asscom2 from any directory on your machine.

```bash
echo "export ASSCOM2_BASE=$(pwd -P)" >> ~/.bashrc # Save installation directory. 
echo "export ASSCOM2_PROFILE=\${ASSCOM2_BASE}/profiles/conda/local" >> ~/.bashrc # Define profile selection. # Different in the development version.
echo "export ASSCOM2_DATABASES=\${ASSCOM2_BASE}/databases" >> ~/.bashrc # Define database base directory.

echo "alias asscom2='conda run --live-stream --name asscom2_launcher snakemake --snakefile \${ASSCOM2_BASE}/snakefile --profile \${ASSCOM2_PROFILE} --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc
source ~/.bashrc

```



## Optionally: Testing the installation

Now you will be able to run asscom2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run asscom2 it will ~~show the message "Pulling singularity image docker://cmkobel/assemblycomparator2." This might take some time depending on your network bandwidth as it downloads a +4GB docker image that contains all the conda environments needed for each analysis.~~ download and install the dependencies of every environment of every rule in the workflow. This might take anywhere from a few minutes to hours. It is generally recommended that you set the channel priority to strict (`conda config --set channel_priority strict`). Unfortunately for one rule "roary" you must disable the channel priority, this is due to Roary not being maintained and some dependencies are not available in the officially specified channels. To disable the channel priority for installing the environment for rule "roary" you can use `conda config --set channel_priority false`. As snakemake installs environments in a random order it is recommended that you first repetitively run `conda config --set channel_priority strict && asscom --until metadata` until all environments but "roary" is installed - Then you can run `conda config --set channel_priority false && asscom --until metadata`. Snakemake installs all environments ahead of time (AOT) regardless of whether you're asking it to run with the "--until" argument.

Nonetheless, if you have any problems installing the conda-based development version, please raise an issue in the githost issues tab.

```bash
# First, enter a dir where some genomes reside.
cd ${ASSCOM2_BASE}/tests/MAGs

# Should take about a minute to complete.
asscom2 --until fast

# Downloads all databases (~ 200 GB).
asscom2 --until downloads

# Run the full pipeline (~ 1 cpu-hour per genome).
asscom2
```