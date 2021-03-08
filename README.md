# assemblycomparator v2

Assemblycomparator is a genomes-to-report pipeline. It is a bit like nullarbor, but it uses genomes (assemblies) instead of reads. 

It is a joint project between Department of Clinical Microbiology Odense at Odense Universityhospital, and Department of Clinical Microbiology Skejby, at Aarhus Universityhospital.


## Installation

 * Install Snakemake: [Official instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) 
 * Set the Snakemake [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) so it matches your cluster setup. 
   * If running locally, no profile is needed
   * If you use slurm, you can utilize the included profile in configs/ pretty much as is.
     * Depending on your cluster configuration, you may also need [drmaa](https://anaconda.org/anaconda/drmaa). 


## What analyses does it do?

As assemblycomparator is currently under development, some analyses are not yes implemented and thouroughly tested.

### For each assembly
  - [ ] [any2fasta](https://github.com/tseemann/any2fasta) (wide input format support)
  - [x] [prokka](https://github.com/tseemann/prokka) (annotation)
  - [ ] [kraken2](https://ccb.jhu.edu/software/kraken2/) (species identification)
  - [x] [mlst](https://github.com/tseemann/mlst) (multi locus sequence typing)
  - [x] [abricate](https://github.com/tseemann/abricate) (virulence/resistance gene identification)
  - [ ] ([Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html)) (Identify possible replication origins, and thereby identify chromids)

  
  
### For each group
  - [x] [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - [ ] [snp-dists](https://github.com/tseemann/snp-dists) (core genome snp-distances)
  - [ ] [panito](https://github.com/sanger-pathogens/panito) (average nucleotide identity
  - [x] [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of core genome)
  - [x] [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - [ ] [IQ-tree](http://www.iqtree.org/) (phylogenetic tree of core genome with bootstrapping)
  - [ ] GC3-profiling ("fingerprinting" of the distribution of GC-content)
  - [ ] Identification of horizontally transferred genes
  - [ ] **A nice report easy to share with your friends and colleagues**
  
  
  
  
