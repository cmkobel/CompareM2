# assemblycomparator2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it takes in genomes (assemblies) instead of reads. Assemblies can come from isolates or metagenomes - as long as they're all prokaryotic.

Assemblycomparator2 works by calling a Snakemake workflow within a conda environment. It performs a palette of analyses on your genomes, and compares them. The main results from these analyses are summarized in a visual html-report that can be easily distributed.

Assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both conda environments and apptainer/singularity/docker images are available for all dependent software to run.


## Usage by examples
Make a directory with the assembly-files you want to investigate with assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
Assemblycomparator2 will then create a sub-directory, named results_ac2/ containing a plethora of analysis results. 

  - Execute a 'dry run'. That is, to show which jobs will run, without submitting any.

    `asscom2 -n`
    
  - Simply, run assemblycomparator on the genomes in the current directory:

    `asscom2`
    
##### A bit more advanced controls 

  - Run analyses that are relevant to metagenomes only:

    `asscom2 --until meta`    
    
  - Execute all jobs until (including) a specific job in the job graph:
    
    `asscom2 --until mlst`
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    `asscom2 --config mlst_scheme=hpylori`
    
  - Select a specific roary blastp-identity: (defaults to 95)

    `asscom2 --config roary_blastp_identity=90`
      


## What analyses does it do?



Below is the graph the shows the dependencies of all possible analyses in assemblycomparator2:

![dag](https://user-images.githubusercontent.com/5913696/236703042-43e1e22c-4013-4c29-a64d-f63fd5b913d5.png)




**Hint:** Use `assemblycomparator2 --until <rulename> [<rulename2>...]` to run specific analyses only. The rulename for each analysis is listed below:

### For each assembly
  - `copy` [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support and validation)
  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) (lengths and GC-content of individual contigs)
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) (generic assembly statistics)
  - `busco` [Busco](https://busco.ezlab.org/) (Estimate assembly completeness and contamination)
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) (Estimate assembly completeness and contamination)
  - `prokka` [prokka](https://github.com/tseemann/prokka) (genomic annotation)
  - `kofam_scan` [KofamScan](https://github.com/takaram/kofam_scan) (annotation of Kegg Orthologous groups, useful for pathway analysis (lacking in report))
  - `dbcan` [dbCAN4](https://github.com/linnabrown/run_dbcan) (annotation of carbohydrate-active "CAZyme" enzymes (lacking in report))
  - `interproscan` [InterProScan](https://github.com/ebi-pf-team/interproscan) (protein function)
  - `abricate` [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
  - `mlst` [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - `kraken2` [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) (Species recognition)
  


### For each group
  - `roary` [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) (core genome pairwise snp-distances (lacking in report))
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of the core genome)
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - **A nice report easy to share with your friends ([demos](https://github.com/cmkobel/assemblycomparator2/tree/master/demo_reports))**


#### Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, kraken2, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, kraken2, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with:
```bash
assemblycomparator2 --until meta
```



## Installation of assemblycomparator2 (asscom2) on Linux

asscom2 can be installed by downloading the code and setting up an alias in your user profile (~/.bashrc) that let's you launch the pipeline from any directory on your machine.

The only requisites for running asscom2 is:
  - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) package manager
  - git distributed version control (can be installed with conda by typing `conda install -c anaconda git`)
  - [apptainer](https://apptainer.org/docs/user/main/quick_start.html#installation-request) container-virtualizer


#### 0) Prerequisites

First, check that you have the prerequisites available on your system:

```bash
which conda && conda --version
which git && git --version
which apptainer && apptainer --version
```

#### 1) Download binary and install launcher environment

Then download the asscom2 pipeline binary and set up an alias in your profile (.bashrc on most linux distributions). Recommended installation directory is in your home directory (\~).

```bash
cd ~
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
        --profile \${ASSCOM2_BASE}/profiles/apptainer/local \
        --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc
source ~/.bashrc

```


## Testing the installation (optional)

Now you will be able to run asscom2. You can use the example data in tests/MAGs to check that everything works. The first time you run asscom2 it will take a long time since it will download a +4GB docker image that contains all the conda packages needed for each analysis.
```bash
cd ${ASSCOM2_BASE}/tests/MAGs
asscom2 --until fast
```

You can also download all the databases to get them ready (~300 GB) with:
```bash
asscom2 --until downloads
```





### Updating an existing installation (optional)

Should you later down the line wish to update the installation, run this command and you should be all set:
```bash

# Pull (download) newest version
cd $ASSCOM2_BASE && git pull

# Install matching version of Snakemake
conda env update --name asscom2_launcher --file environment.yaml


```





## Future functionality 

In the future we might add some of the following pieces of software into assemblycomparator2.

**Sample basis**

  - [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) (Identify possible replication origins, and thereby help identify chromids)
  - [RFplasmid](https://github.com/aldertzomer/RFPlasmid) (Identify plasmids using the pentamer-random-forest method)
  - [Kaptive](https://github.com/katholt/Kaptive) (surface polysaccharide loci for Klebsiella and Acinetobacter baumannii) 
  - [mash screen](https://mash.readthedocs.io/en/latest/tutorials.html) (recognition of plasmids-of-interest)


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
