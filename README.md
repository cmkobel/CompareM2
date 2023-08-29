# assemblycomparator2

[![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/assemblycomparator2?label=docker%20pulls%202023)](#installation-of-assemblycomparator2-on-linux)

assemblycomparator2 (asscom2) is a genomes-to-report pipeline. It is a bit like nullarbor, but it takes in genomes (assemblies) instead of reads. Assemblies can come from isolates or metagenomes - as long as they're all prokaryotic.

assemblycomparator2 works by calling a Snakemake workflow within a conda environment. It performs a palette of 16 analyses on your genomes, and compares them. The main results from these analyses are summarized in a visual portable .html-document report that can be easily shared. This report is generated even if a few jobs in the pipeline fail.

assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  apptainer/singularity/docker images and conda environment definitions are available for all dependent software to run.


## Usage examples

Make a directory with the assembly-files you want to investigate with assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
assemblycomparator2 will then create a sub-directory, named results_ac2/ containing a plethora of analysis results. 
  
  - Execute a 'dry run'. That is, to show what will be run without actually doing it.

    ```bash
    asscom2 --dry-run
    ```

  - Run assemblycomparator2 on the genomes in the current directory:

    ```bash
    asscom2
    ```
    

##### A bit more advanced controls 

  - Run analyses that are relevant to metagenomes only:

    ```bash
    asscom2 --until meta
    ```
    
  - Execute all jobs until (including) a specific job in the job graph:
    
    ```bash
    asscom2 --until mlst
    ```
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    ```bash
    asscom2 --config mlst_scheme=hpylori
    ```
    
  - Select a specific roary blastp-identity: (default is 95)

    ```bash
    asscom2 --config roary_blastp_identity=90
    ```
      


## What analyses does it do?

Below is the graph the shows the order of execution of all possible analyses in assemblycomparator2:

![dag](https://github.com/cmkobel/assemblycomparator2/assets/5913696/5edc09e0-fe5f-4d7f-aeff-b72acd728ba7)



**Hint:** Use `asscom2 --until <rulename> [<rulename2>...]` to run one or several specific analyses only. The rulename for each analysis to pick is listed below:

### For each assembly
  - `copy` [any2fasta](https://github.com/tseemann/any2fasta) Wide input format support and validation.
  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) Lengths and GC-content of individual contigs.
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) Generic assembly statistics.
  - `busco` [Busco](https://busco.ezlab.org/) Estimate assembly completeness and contamination.
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) Estimate assembly completeness and contamination.
  - `prokka` [prokka](https://github.com/tseemann/prokka) Genomic annotation.
  - `kofam_scan` [KofamScan](https://github.com/takaram/kofam_scan) Annotation of Kegg Orthologous groups.
  - `diamond_kegg` [diamond](https://github.com/bbuchfink/diamond) Run prokka-called proteins through the checkm2 database (uniref100 with KEGG-KOs).
  - `kegg_pathway` [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) KEGG ortholog-based hypergeometric pathway enrichment analysis.
  - `dbcan` [dbCAN4](https://github.com/linnabrown/run_dbcan) Annotation of carbohydrate-active "CAZyme" enzymes (lacking in report).
  - `interproscan` [InterProScan](https://github.com/ebi-pf-team/interproscan) Protein function using Tigrfam, Hamap and Pfam.
  - `abricate` [abricate](https://github.com/tseemann/abricate) Virulence and resistance gene identification.
  - `mlst` [mlst](https://github.com/tseemann/mlst) Multi locus sequence typing.
  - `kraken2` [kraken2](https://ccb.jhu.edu/software/kraken2/) Species identification.
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) Species recognition.
  

### For each group
  - `roary` [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) (core genome pairwise snp-distances (lacking in report))
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of the core genome)
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - **A nice report easy to share with your friends ([demos](https://github.com/cmkobel/assemblycomparator2/blob/master/readme-demos.md))**


#### Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, kraken2, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, kraken2, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with:
```bash
asscom --until meta

asscom2 --until isolate
```



## Installation of assemblycomparator2 on Linux

assemblycomparator2 can be installed by downloading the code and setting up an alias in your user profile that let's you launch the pipeline from any directory on your machine.

The only requisites for running assemblycomparator2 is:
  - [*conda*](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) package manager
  - *git* distributed version control (can be installed with conda by typing `conda install -c anaconda git`)
  - [*apptainer*](https://apptainer.org/docs/user/main/quick_start.html#installation-request) container-virtualizer


#### 0) Prerequisites

First, check that you have the prerequisites available on your system:

```bash
which conda && conda --version
which git && git --version
which apptainer && apptainer --version
```

#### 1) Download pipeline and set up the launcher environment

Then download the assemblycomparator2 pipeline and set up an alias in your profile (.bashrc on most linux distributions). Proposed installation directory is in your home directory (\~).

```bash
cd ~
git clone https://github.com/cmkobel/assemblycomparator2.git asscom2
cd asscom2
conda env create --name asscom2_launcher --file environment.yaml # Installs snakemake and mamba in an environment named "asscom2_launcher".

```


#### 2) Alias

Finally, define the alias that will be used to launch asscom2 from any directory on your machine.

```bash
echo "export ASSCOM2_BASE=$(pwd -P)" >> ~/.bashrc # Save installation directory. 
echo "export ASSCOM2_PROFILE=\${ASSCOM2_BASE}/profiles/apptainer/local" >> ~/.bashrc # Save profile selection.
echo "alias asscom2='conda run --live-stream --name asscom2_launcher snakemake --snakefile \${ASSCOM2_BASE}/snakefile --profile \${ASSCOM2_PROFILE} --configfile \${ASSCOM2_BASE}/config.yaml'" >> ~/.bashrc
source ~/.bashrc

```


## Testing the installation (optional)

Now you will be able to run asscom2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run asscom2 it will show the message "Pulling singularity image docker://cmkobel/assemblycomparator2." This might take some time depending on your network bandwidth as it downloads a +4GB docker image that contains all the conda environments needed for each analysis.
```bash
# First, enter a dir where some genomes reside.
cd ${ASSCOM2_BASE}/tests/MAGs

# Should take about a minute to complete.
asscom2 --until fast

# Downloads all databases (~ 200 GB).
asscom2 --until downloads

# Run the full pipeline (~ 1 cpu-hour per genome).
asscom2
```





### Updating an older installation (optional)

If you want to make sure that you're running the latest version of assemblycomparator2, you can run these commands to update it:
```bash

# Pull (download) newest version
cd $ASSCOM2_BASE && git pull

# Install matching version of Snakemake
conda env update --name asscom2_launcher --file environment.yaml


```





## Future functionality 

In the future we might add some of the following software packages into assemblycomparator2.

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

If you use assemblycomparator2, you can support further funding by bumping the citation count on this one:

  - Kobel CM. *assemblycomparator2* **GitHub** https://github.com/cmkobel/assemblycomparator2

assemblycomparator2 would not have existed, if it hadn't been for the integrated software packages and their databases.

  - W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962. 
  - Richard Challis. (2017). rjchallis/assembly-stats 17.02 (17.02). Zenodo. https://doi.org/10.5281/zenodo.322347
  - Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
  - Chklovski, A., Parks, D.H., Woodcroft, B.J. et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods 20, 1203–1212 (2023). https://doi.org/10.1038/s41592-023-01940-w
  - Seemann T., Prokka: rapid prokaryotic genome annotation, Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063
  - Takuya Aramaki and others, KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold, Bioinformatics, Volume 36, Issue 7, April 2020, Pages 2251–2252, https://doi.org/10.1093/bioinformatics/btz859
  - Yin Y*, Mao X*, Yang JC, Chen X, Mao F and Xu Y, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Res. 2012 Jul;40(Web Server issue):W445-51 
  - P. Jones et al., Bioinformatics (2014), PMID: 24451626
  - Seemann T, Abricate, Github https://github.com/tseemann/abricate
  - Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018 Sep 24;3:124. doi: 10.12688/wellcomeopenres.14826.1. PMID: 30345391; PMCID: PMC6192448.
  - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
  - Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics. 2019 Nov 15;36(6):1925–7. doi: 10.1093/bioinformatics/btz848. Epub ahead of print. PMID: 31730192; PMCID: PMC7703759.
