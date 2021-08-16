# assemblycomparator2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it takes in genomes (assemblies) instead of reads. 

It is a joint project between Department of Clinical Microbiology Odense at Odense Universityhospital, and Department of Clinical Microbiology Skejby, at Aarhus Universityhospital.

## Usage
Make a directory with the assembly-files you want to investigate with assemblycomparator2. 
Go into that directory in the terminal, and run the command `assemblycomparator2_slurm` or `assemblycomparator2_local`. 
assemblycomparator2 will then create a sub-directory containing a plethora of analyses. 

#### Some useful commands
  - Execute a 'dry run'. That is, show the jobs which will run, without triggering the computation:
    
    `assemblycomparator2_slurm -n`
    
  - Simply, run assemblycomparator on the genomes in the current directory:

    `assemblycomparator2_slurm`
    
  - If you're not sure your internet connection to the cluster will last for the full assemblycomparator2 run, put a `&` in the end.
  
    `assemblycomparator2_slurm &`
    
  - Execute all jobs up until (inclusive of) a specific job in the job graph:
    
    `assemblycomparator2_slurm --until mlst`
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    `assemblycomparator2_slurm --config mlst_scheme=hpylori`
    
  - Select a specific roary blastp-identity: (defaults to 95)

    `assemblycomparator2_slurm --config roary_blastp_identity=90`
    
  - Rerun a specific rule, (might be necessary if some parts of the report is missing):

    `assemblycomparator2_slurm -R report`
    
    





## What analyses does it do?

### For each assembly
  - [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support)
  - [prokka](https://github.com/tseemann/prokka) (annotation)
  - [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
  - [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) (generic assembly statistics)


### For each group
  - [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - [snp-dists](https://github.com/tseemann/snp-dists) (core genome pairwise snp-distances)
  - [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of the core genome)
  - [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - **A nice report easy to share with your friends ([demo](https://github.com/cmkobel/assemblycomparator2/raw/master/tests/E._faecium/report_E._faecium.html.zip))**


Below is a snakemake exported directed graph of the rules involved:
![dag](https://user-images.githubusercontent.com/5913696/123612417-cd402780-d802-11eb-81e8-617356fb62cf.png)






## Installation

Assemblycomparator2 needs Snakemake and the dependencies which can be needed for running on your specific setup. I.e. DRMAA for Slurm-mananged HPC's.
You can either follow the [official Snakemake instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or use our guide below:
* Decide where you want to install assemblycomparator2
   ``` 
   ASSCOM2_BASE=~/assemblycomparator2
   
   mkdir -p $ASSCOM2_BASE
    
   # And save it into your .bashrc
   echo "export ASSCOM2_BASE=$ASSCOM2_BASE" >> ~/.bashrc 
    
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
 * You have to decide whether you want to use Singularity (recommended if possible) or Conda for package management.

   
   
### For HPC's with Slurm using Singularity
   ```
   # Main alias for running assemblycomparator2
   echo "alias assemblycomparator2_slurm='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/slurm/ --cluster-config ${ASSCOM2_BASE}/configs/slurm/slurm.yaml --use-singularity  --singularity-prefix ${ASSCOM2_BASE}/singularity_images --configfile ${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc
    
   ```
   
   
### For local setups using Conda
   ```
   # Main alias for running assemblycomparator2
   echo "alias assemblycomparator2_local='conda activate assemblycomparator2; snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/local/ --use-conda --configfile ${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc
   
   # Set the SNAKEMAKE_CONDA_PREFIX-variable, so the package installations can be reused between runs.
   echo "export SNAKEMAKE_CONDA_PREFIX=${ASSCOM2_BASE}/conda_base" >> ~/.bashrc 
    
   ```
   
   
 * assemblycomparator2 supports Kraken2. If you already have a local copy of a kraken database, you can set the `ASSCOM2_KRAKEN_DB` system variable to its path. If you don't have a local copy, assemblycomparator2 comes handy with some scripts for setting up Kraken2 and Mashscreen. There are two scripts for Kraken2; one small "Standard" (8GB) and one huge "PlusPF" (50GB).
   ```
   # Kraken2
   $ASSCOM2_BASE/scripts/set-up-kraken2_Standard_8GB.sh
   # or:
   $ASSCOM2_BASE/scripts/set-up-kraken2_PlusPF_50GB.sh
   
   # Mashscreen
   $ASSCOM2_BASE/scripts/set_up_mashscreen.sh
    
   ```
 * When you have completed all installation steps, you should read the settings into global system memory. After this, you are ready to test the installation.
   ```
   source ~/.bashrc
    
   ```
   
### Testing installation

assemblycomparator2 comes with a handful of E. faecium assemblies (illumina/skesa) which can be used to check that everything works as expected. In order to run this test, simply go into the location of these assemblies, and run the `assemblycomparator2`-command
   ```
   cd ${ASSCOM2_BASE}/tests/E._faecium_plasmids
   assemblycomparator2
    
   ```

If you encounter problems testing your installation, please refer to the issues tab of this repository.

   
   
### Updating an existing installation

Simply run this command, and you should be all set:
```
cd $ASSCOM2_BASE && git pull

# You might also want to update snakemake
conda env update --name assemblycomparator2 --file environment.yaml
 
```
Note: If new databases have been added to kraken or mashscreen, you can rerun the above-mentioned set_up_*.sh-scripts.




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
  - [panito](https://github.com/sanger-pathogens/panito) (average nucleotide identity)
  - [GenAPI](https://github.com/MigleSur/GenAPI) (alternative to roary)


  
  
  
