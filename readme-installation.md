
--------------------
## Installation

Dependencies:
  - Snakemake 
  - Either Apptainer or Conda

Assemblycomparator2 needs Snakemake to run. By default, assemblycomparator2 uses singularity to run containerized environments (images) that are automatically pulled from dockerhub. Alternatively, it is also possible to automatically install and use corresponding conda environments, but this is recommended only for developers).

### 1. Set up the assemblycomparator2 base directory

Firstly, define the base directory for assemblycomparator2 in $ASSCOM2_BASE. If neccessary, you can change it to anything you'd like.

```bash
# Define base path.
ASSCOM2_BASE=~/assemblycomparator2

# Create base directory.
mkdir -p $ASSCOM2_BASE

# Save it into your ~/.bashrc
echo "export ASSCOM2_BASE=$ASSCOM2_BASE" >> ~/.bashrc  
```

### 2. Download (clone) assemblycomparator2

Clone the assemblycomparator2 git-repository into that base. If you don't have git, you can also just download the files.

```bash
git clone https://github.com/cmkobel/assemblycomparator2.git $ASSCOM2_BASE

# Developers can preferably use the git protocol:
# git clone git@github.com:cmkobel/assemblycomparator2.git $ASSCOM2_BASE
```

### 3. Installing Snakemake

The recommended way of installing Snakemake is to make a conda environment that contains the exact Snakemake version needed. This environment is only used to start the interactive snakemake job scheduler.

```bash

#  It is recommended to use the strict channel policy in conda. 
conda config --set channel_priority strict

# It is also recommended that you use mamba to install any package within any conda environment
conda install -c conda-forge mamba

# Install Snakemake in a environment by using the bundled yaml-file.
cd $ASSCOM2_BASE && mamba env create -n assemblycomparator2 -f environment.yaml

```

You don't need to do any manual activation of the conda environment. As assemblycomparator2 is wrapped within conda run, it automatically activates the correct environment for each analysis.
   
### 4. Set up the alias 

#### HPC workload manager or local execution?
<table><tr><td>
If you don't know what a HPC workload manager is, or don't want to use it, skip directly to <a href="https://github.com/cmkobel/assemblycomparator2#41-setting-up-assemblycomparator2-to-use-apptainer-on-a-local-workstation">subheading 4.1</a>.
</td></tr></table>

Assemblycomparator2 supports both Apptainer (formerly known as Singularity) and Conda. When using Apptainer, the images are automatically pulled from Dockerhub. When using Conda, the environments are automatically solved and installed into your ~/.asscom2/conda_base/ directory. <u>It is strongly recommended that you use Apptainer</u>, as Conda can be very slow at solving environments on HPCs with congested network drives, even when using mamba.

Assemblycomparator2 supports both running directly on the "local" machine (i.e. workstation), or through a HPC workload manager such as "SLURM", "PBS" or others. 

Your decision is to pick the profile within the profiles/ directory that describes the way that you want to run assemblycomparator2. Below is the file tree that shows the options that are available. Please note that you should give the path to the folder containing the config.yaml file that you want to use. If there isn't a profile that matches your system you should feel free to create one by using one of the bundled as template.

```bash
#> profiles/
#> ├── apptainer/
#> │   ├── local/
#> │   ├── pbs-qut-lyra/
#> │   ├── slurm-au-genomedk/
#> │   ├── slurm-nmbu-orion/
#> │   └── slurm-sigma2-saga/
#> └── conda/
#>     ├── local/
#>     ├── pbs-qut-lyra/
#>     ├── slurm-au-genomedk/
#>     ├── slurm-nmbu-orion/
#>     └── slurm-sigma2-saga/
```


Save this decision in the ASSCOM2_PROFILE variable. Below is an example shown of setting assemblycomparator2 up to use apptainer on a local workstation, which is the recommended option for new users. 

#### 4.1 Setting up assemblycomparator2 to use apptainer on a local workstation.

```bash
# Define the Snakemake profile to use
ASSCOM2_PROFILE=${ASSCOM2_BASE}/profiles/apptainer/local

# Save it into your ~/.bashrc
echo "export ASSCOM2_PROFILE=$ASSCOM2_PROFILE" >> ~/.bashrc 
```

Finally, after defining the ASSCOM2_PROFILE, we can define the alias that you will use to run assemblycomparator2.

#### 4.2 Finalizing install by defining the alias that starts assemblycomparator2
```bash
# Main alias for running assemblycomparator2
echo "alias assemblycomparator2='conda run --live-stream --name assemblycomparator2 \
    snakemake \
        --snakefile \${ASSCOM2_BASE}/snakefile \
        --profile \${ASSCOM2_PROFILE} \
        --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc


# Then source your ~/.bashrc
source ~/.bashrc

```

Now, the installation of assemblycomparator2 is complete.
   
### Testing the assemblycomparator2 installation (optional)

Assemblycomparator2 comes with a few directories with test assemblies which can be used to check that everything works as expected. In order to run this test, simply go into the location of these assemblies, and run the `assemblycomparator2`-command
```bash
cd ${ASSCOM2_BASE}/tests/E._faecium_plasmids

assemblycomparator2 
```

If you encounter problems installing, testing or using assemblycomparator2, please an issue on the [github issue tracker](https://github.com/cmkobel/assemblycomparator2/issues). No issues are too small.

   
   