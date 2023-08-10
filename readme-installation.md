## Installation of assemblycomparator2 (asscom2) on Linux

asscom2 can be installed by downloading the code and setting up an alias in your user profile (~/.bashrc) that let's you launch the pipeline from any directory on your machine.

The only requisites for running asscom2 is:
  - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) package manager
  - git (can be installed with conda by typing `conda install -c anaconda git`) distributed version control
  - [apptainer](https://apptainer.org/docs/user/main/quick_start.html#installation-request) container-virtualizer


#### 0) Prerequisites

First, check that you have the prerequisites available on your system:

```bash
which conda && conda --version
which git && git --version
which apptainer && apptainer --version
```

#### 1) Download binary and install launcher environment

Then download the asscom2 pipeline binary and set up an alias in your profile (\~/.bashrc on most linux distributions). Recommended installation directory is your home directory (\~).

```bash
# cd ~
git clone https://github.com/cmkobel/assemblycomparator2.git asscom2
cd asscom2
conda env create --name asscom2_launcher --file environment.yaml # Installs snakemake and mamba in an environment named "asscom2_launcher".

```


#### 2) Alias

Finally, define the alias that will be used to launch asscom2 from any directory on your machine.

```bash
echo "export ASSCOM2_BASE=$(pwd)" >> ~/.bashrc
echo "alias asscom2='conda run --live-stream --name asscom2_launcher \
    snakemake \
        --snakefile \${ASSCOM2_BASE}/snakefile \
        --profile \${ASSCOM2_BASE}/apptainer/local \
        --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc

```







