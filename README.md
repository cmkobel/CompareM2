# assemblycomparator v2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it uses genomes (assemblies) instead of reads. 

It is a joint project between Department of Clinical Microbiology Odense at Odense Universityhospital, and Department of Clinical Microbiology Skejby, at Aarhus Universityhospital.


## Installation

Assemblycomparator2 needs Snakemake and the dependencies which can be needed for running on your specific setup. I.e. DRMAA for Slurm-mananged HPC's.
You can either follow the [official Snakemake instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or use our guide below:
* Decide where you want to install assemblycomparator2
   ``` 
   ASSCOM2_BASE=~/assemblycomparator2
   
   # And save it into your .bashrc
   echo "export ASSCOM2_BASE=$ASSCOM2_BASE" >> ~/.bashrc && source ~/.bashrc
   ```
 * Clone the assemblycomparator2 GitHub-repository into that base
   ```
   git clone https://github.com/cmkobel/assemblycomparator2.git $ASSCOM2_BASE
   
   # Hint: If you haven't already installed Snakemake and its dependencies, you can do it easily now:
   conda env create -f environment.yml 
   ```
   
 * Set an alias that makes it easy to run assemblycomparator2 from anywhere in your filesystem
   ```
   # For local setups:
   echo "alias assemblycomparator2='snakemake --snakefile ${ASSCOM2_BASE}/snakefile'" >> ~/.bashrc
   
   # For HPC's with Slurm:
   echo "alias assemblycomparator2='snakemake --snakefile ${ASSCOM2_BASE}/snakefile --profile ${ASSCOM2_BASE}/configs/slurm/ --cluster-config ${ASSCOM2_BASE}/configs/cluster.yaml'" >> ~/.bashrc
   ```
   
 * Optionally: Consider running the Kraken2 and mash screen set up scripts:
   ```
   $ASSCOM2_BASE/scripts/set_up_kraken2.sh
   $ASSCOM2_BASE/scripts/set_up_mashscreen.sh
   ```
   
   
   
   
### Updating an existing installation

Simply run this command, and you should be all set:
```
cd $ASSCOM2_BASE
git pull
```
Note: If new database have been added to kraken or mashscreen, you can rerun the above-mentioned set_up_*.sh-scripts.



## What analyses does it do?

As assemblycomparator is currently under development, some analyses are not yet implemented and thoroughly tested.

### For each assembly
  - [ ] [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support)
  - [x] [prokka](https://github.com/tseemann/prokka) (annotation)
  - [ ] [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - [x] [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - [x] [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
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
  
  
  
  
