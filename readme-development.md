# Development information

This file describes how to install assemblycomparator2 for development purposes. Since you will be defining and configuring a conda environment that will eventually be added to the apptainer-compatible dockerfile, you must run everything through conda. Once you make a pull request and it is accepted, your conda environment will be added to the official docker image which is then published with the next version number.

Assemblycomparator2 only works on Linux. There is no plan to support any other operating system.

## Installation of assemblycomparator2

assemblycomparator2 can be installed by downloading the code and setting up an alias in your user profile that let's you launch the pipeline from any directory on your machine.

The only requisites for running assemblycomparator2 is:
  - conda/[mamba](https://github.com/conda-forge/miniforge#install) package manager
  - *git* distributed version control

#### 0) Prerequisites

First, check that you have the prerequisites available on your system. 

```bash
which mamba && mamba --version
which git && git --version
```

#### 1) Download pipeline and set up the launcher environment

Then download the assemblycomparator2 pipeline and set up an alias in your profile (.bashrc on most linux distributions). Proposed installation directory is in your home directory (\~).

```bash
cd ~ # Enter the directory where you want to install assemblycomparator2.
git clone https://github.com/cmkobel/assemblycomparator2.git asscom2
cd asscom2
mamba env create -y -f environment.yaml -n asscom2_dev

```

#### 2) Force using Conda to create and run individual tools.

Using the profile on this path disables the Docker compatible apptainer image.

_Skippable_: If you are planning to make changes in Assemblycomparator2 that do not add or modify individual Conda environments, you can skip this step.

```bash
export ASSCOM2_PROFILE="${ASSCOM2_BASE}/profile/conda/default"

```


#### 3) Activate conda env and run

Finally, you can run the pipeline with the following code:
```bash
conda activate asscom2_dev
unzip tests/E._faecium/fna.zip # Or gather your own relevant files for testing.
./asscom2 --config input_genomes="*.fna" --until fast

```

On the first run you should see Snakemake creating the conda environments.

#### 4) Done

Happy hacking!
If in doubt, please make an issue to check if there is interest in developing a new function.
