# assemblycomparator2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it takes in genomes (assemblies) instead of reads. 

It is a joint project between Department of Clinical Microbiology Odense at Odense Universityhospital, and Department of Clinical Microbiology Skejby, at Aarhus Universityhospital.

## Usage
Make a directory with a few fasta files you want to investigate and possibly compare. 
Go into that directory in the terminal, and run the command `assemblycomparator2`. 
assemblycomparator2 will then create a sub-directory containing a plethora of analyses. 

We are also working on finishing an automatically generated report, which will make it easy to explore these analyses without rummaging around on the filesystem.
That is it!



## Installation

Assemblycomparator2 needs Snakemake and the dependencies which can be needed for running on your specific setup. I.e. DRMAA for Slurm-mananged HPC's.
You can either follow the [official Snakemake instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or use our guide below:
* Decide where you want to install assemblycomparator2
   ``` 
   ASSCOM2_BASE=~/assemblycomparator2
   
   # And save it into your .bashrc
   echo "export ASSCOM2_BASE=$ASSCOM2_BASE" >> ~/.bashrc 
   
   # If you are going to use conda, it is a good idea to set the SNAKEMAKE_CONDA_PREFIX-variable, so the package installations can be reused between runs.
   echo "export SNAKEMAKE_CONDA_PREFIX=${ASSCOM2_BASE}/conda_base" >> ~/.bashrc 
   ```
 * Clone the assemblycomparator2 GitHub-repository into that base
   ```
   git clone https://github.com/cmkobel/assemblycomparator2.git $ASSCOM2_BASE
   
   # Optionally use the git protocol:
   # git clone git@github.com:cmkobel/assemblycomparator2.git $ASSCOM2_BASE
   
   # Hint: If you haven't already installed Snakemake and its dependencies, you can do it easily now. (Might take a few minutes):
   cd $ASSCOM2_BASE && conda env create -f environment.yaml 
   ```
   
 * Set an alias that makes it easy to run assemblycomparator2 from anywhere in your filesystem
   ```
   # For local setups (using Conda for jobs):
   echo "alias assemblycomparator2='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/local/ --use-conda'" >> ~/.bashrc
   
   # For HPC's with Slurm (using Singularity for jobs):
   echo "alias assemblycomparator2='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/slurm/ --cluster-config ${ASSCOM2_BASE}/configs/cluster.yaml --use-singularity'" >> ~/.bashrc
   ```
   Hint: You can interchange `--use-conda` and `--use-singularity` for changing how assemblycomparator2 runs the jobs. Please note that running assemblycomparator2 locally with conda is not fully developed or tested, and has a high probability of failing. If you have access to Singularity, use it.
   
   * assemblycomparator2 supports Kraken2. If you already have a local copy of a kraken database, you can set the `ASSCOM2_KRAKEN_DB` system variable to its path. If you don't have a local copy, assemblycomparator2 comes handy with some scripts for setting up Kraken2 and Mashscreen. There are two scripts for Kraken2; one small "Standard" (8GB) and one huge "PlusPF" (50GB).
   ```
   # Kraken2
   $ASSCOM2_BASE/scripts/set-up-kraken2_Standard_8GB.sh
   # or:
   $ASSCOM2_BASE/scripts/set-up-kraken2_PlusPF_50GB.sh
   
   # Mashscreen
   $ASSCOM2_BASE/scripts/set_up_mashscreen.sh
   ```
   * When you have completed all installation steps, you should read the settings into global system memory.
   ```
   source ~/.bashrc
   ```
   
   
   
### Updating an existing installation

Simply run this command, and you should be all set:
```
cd $ASSCOM2_BASE && git pull

# You might also want to update snakemake
conda env update --name assemblycomparator2 --file environment.yaml
```
Note: If new databases have been added to kraken or mashscreen, you can rerun the above-mentioned set_up_*.sh-scripts.



## What analyses does it do?

As assemblycomparator is currently under development, some analyses are not yet implemented or thoroughly tested.

### For each assembly
  - [x] [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support)
  - [x] [prokka](https://github.com/tseemann/prokka) (annotation)
  - [x] [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - [x] [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - [x] [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
  - [x] [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) (generic assembly statistics)
  - [ ] [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) (Identify possible replication origins, and thereby identify chromids)
  - [ ] [RFplasmid](https://github.com/aldertzomer/RFPlasmid) (Identify plasmids using the pentamer-random-forest method)
  - [ ] [Kaptive](https://github.com/katholt/Kaptive) (surface polysaccharide loci for Klebsiella and Acinetobacter baumannii)

  
  
### For each group
  - [x] [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - [ ] [GenAPI](https://github.com/MigleSur/GenAPI) (alternative to roary)
  - [ ] [snp-dists](https://github.com/tseemann/snp-dists) (core genome snp-distances)
  - [ ] [panito](https://github.com/sanger-pathogens/panito) (average nucleotide identity
  - [x] [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of core genome)
  - [x] [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - [ ] [IQ-tree](http://www.iqtree.org/) (phylogenetic tree of core genome with bootstrapping)
  - [ ] GC3-profiling ("fingerprinting" of the distribution of GC-content)
  - [ ] Identification of horizontally transferred genes
  - [ ] **A nice report easy to share with your friends and colleagues**
  
  
  
  
