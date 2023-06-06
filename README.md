# assemblycomparator2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it takes in genomes (assemblies) instead of reads. Assemblies can come from isolates or metagenomes - as long as they're all prokaryotic.

Assemblycomparator2 works by calling a Snakemake workflow within a conda environment. It performs a palette of analyses on your genomes, and compares them. The main results from these analyses are summarized in a visual html-report that can be easily distributed.

Assemblycomparator2 can be run either on a local workstation (>64GiB RAM), or a HPC (high performance computing) cluster. Both conda environments and apptainer/singularity/docker images are available for all dependent software to run.


## Usage by examples
Make a directory with the assembly-files you want to investigate with assemblycomparator2. 
Go into that directory in the terminal, and run the command `assemblycomparator2`. 
Assemblycomparator2 will then create a sub-directory, named results_ac2/ containing a plethora of analysis results. 

  - Execute a 'dry run'. That is, to show which jobs will run, without submitting any.

    `assemblycomparator2 -n`
    
  - Simply, run assemblycomparator on the genomes in the current directory:

    `assemblycomparator2`
    
##### A bit more advanced controls 
    
  - Execute all jobs until (including) a specific job in the job graph:
    
    `assemblycomparator2 --until mlst`
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    `assemblycomparator2 --config mlst_scheme=hpylori`
    
  - Select a specific roary blastp-identity: (defaults to 95)

    `assemblycomparator2 --config roary_blastp_identity=90`
      


## What analyses does it do?



Below is the graph the shows the dependencies of all possible analyses in assemblycomparator2:

![dag](https://user-images.githubusercontent.com/5913696/236703042-43e1e22c-4013-4c29-a64d-f63fd5b913d5.png)




**Hint:** Use `assemblycomparator2 --until <rulename> [<rulename2>...]` to run specific analyses only. The rulename for each analysis is listed below:

### For each assembly
  - `copy` [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support)
  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) (lengths and GC-content of individual contigs)
  - `prokka` [prokka](https://github.com/tseemann/prokka) (annotation)
  - `kraken2` [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - `mlst` [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - `abricate` [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) (generic assembly statistics)
  - \#[clusterProfiler KEGG](https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html) (pathway enrichment analysis)
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) (Species recognition)
  - `busco` [Busco](https://busco.ezlab.org/) (Estimate assembly completeness and contamination)
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) (Estimate assembly completeness and contamination)

  


### For each group
  - `roary` [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) (core genome pairwise snp-distances)
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of the core genome)
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - **A nice report easy to share with your friends ([demos](https://github.com/cmkobel/assemblycomparator2/tree/master/demo_reports))**


#### Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, kraken2, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, kraken2, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with:
```bash
assemblycomparator2 --until meta
```



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

   
   
### Updating an existing installation (optional)

If you should -later down the line- wish to update the installation, run this command and you should be all set:
```bash

# Pull (download) newest version
cd $ASSCOM2_BASE && git pull

# Install matching version of Snakemake
conda env update --name assemblycomparator2 --file environment.yaml


```





## Future functionality 

In the future we might add some of the following pieces of software into assemblycomparator2.

**Sample basis**

  - [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) (Identify possible replication origins, and thereby help identify chromids)
  - [RFplasmid](https://github.com/aldertzomer/RFPlasmid) (Identify plasmids using the pentamer-random-forest method)
  - [Kaptive](https://github.com/katholt/Kaptive) (surface polysaccharide loci for Klebsiella and Acinetobacter baumannii) 
  - [mash screen](https://mash.readthedocs.io/en/latest/tutorials.html) (recognition of plasmids-of-interest)
  - [SingleM](https://wwood.github.io/singlem/)


**Batch basis**

  - [IQ-tree](http://www.iqtree.org/) (phylogenetic tree of core genome with bootstrapping)
  - GC3-profiling ("fingerprinting" of the distribution of GC-content)
  - Identification of horizontally transferred genes?
  - [GenAPI](https://github.com/MigleSur/GenAPI) (alternative to roary)


---

Development is active and will continue.

Assemblycomparator2 genomes to report pipeline. Copyright (C) 2019-2023 [Carl M. Kobel](https://github.com/cmkobel) [GNU GPL v3](https://github.com/cmkobel/assemblycomparator2/blob/master/LICENSE)
  
  

## Citation

Kobel CM. *assemblycomparator2* **GitHub** https://github.com/cmkobel/assemblycomparator2
