# Development information

This file describes how to install comparem2 for development purposes. If you want to add new analysis rules you will be defining and configuring a conda environment that will eventually be added to the apptainer-compatible dockerfile, you must run everything through Conda. Once you make a pull request and it is accepted, your conda environment will be added to the official docker image which is then published with the next version number.

CompareM2 only works on Linux.

## Installation of comparem2

comparem2 can be installed by downloading the code and setting up an alias in your user profile that let's you launch the pipeline from any directory on your machine.

The only requisites for running comparem2 is:
  - conda compatible package manager (We recommend [mamba from miniforge](https://github.com/conda-forge/miniforge#install))
  - *git* distributed version control

#### 0) Prerequisites

First, check that you have the prerequisites available on your system. 

```bash
which mamba && mamba --version
which git && git --version
```

#### 1) Download pipeline and set up the launcher environment

Then download the comparem2 pipeline and set up an alias in your profile (.bashrc on most linux distributions). Proposed installation directory is in your home directory (\~).

```bash
cd ~ # Enter the directory where you want to install comparem2.
git clone https://github.com/cmkobel/comparem2.git comparem2
cd comparem2
mamba env create -y -f environment.yaml -n comparem2_dev

```

If you are planning to make a pull request, you can also clone your personal fork instead of cloning the official repository.

#### 2) Force using Conda to create and run individual tools.

If you have Apptainer on your machine, you should use this specific profile to force using Conda instead.

_Skippable_: If you are planning to make changes in CompareM2 that do not add or modify individual Conda environments, you can skip this step.

```bash
export COMPAREM2_PROFILE="$(realpath profile/conda/default)"

```

`realpath` makes sure that you get the absolute path. This is helpful, should you later change directory.

You can read more about CompareM2 environment variables here: https://comparem2.readthedocs.io/en/latest/10%20installation/#advanced-configuration

#### 3) Activate conda env and run

Finally, you can run the pipeline with the following code:
```bash
conda activate comparem2_dev
unzip tests/E._faecium/fna.zip # Or gather your own relevant files for testing.

./comparem2 --config input_genomes="*.fna" --until fast --dry-run

```

Issuing a _dry run_ is a good way to get started, checking that everything is OK.

#### 4) Done

Happy hacking!

