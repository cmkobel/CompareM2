# assemblycomparator2

[![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/assemblycomparator2?label=docker%20pulls%202023)](#installation-of-assemblycomparator2-on-linux) [![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/cmkobel/assemblycomparator2)](https://hub.docker.com/r/cmkobel/assemblycomparator2) [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/assemblycomparator2?label=Bioconda%20downloads&color=%2300CC00)](https://anaconda.org/bioconda/assemblycomparator2) [![conda build](https://img.shields.io/conda/v/bioconda/assemblycomparator2)](https://anaconda.org/bioconda/assemblycomparator2)

Assemblycomparator2 (asscom2) is a genomes-to-report pipeline. It accepts prokaryotic genomic assemblies and compares them in many different ways. Designed to analyze both isolates and metagenomes.

![2ezgif-3-6b2c1d8231](https://github.com/cmkobel/assemblycomparator2/assets/5913696/623f6b42-2de6-457c-8f0d-3b3e5d646967)



Assemblycomparator2 integrates several analyses that yield scientific results about genomic assemblies on several levels: Quality control, annotation, function and species calling as well as comparative analyses like computation of core/pan genomes and phylogenetics. 

Assemblycomparator2 works by calling a Snakemake workflow that can be easily modified to use different parameters for the  underlying tools.

All results are dynamically integrated in a compact portable report .html-document that emphasizes the central results and can be easily shared. This report is generated even if a few jobs in the pipeline fail.

Assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  apptainer/singularity/docker images and conda environment definitions are available for all dependent software to run.

Assemblycomparator2 will soon be renamed to "Proknome".

## Usage examples

Make a directory with the assembly-files you want to investigate with assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
assemblycomparator2 will then create a sub-directory, named "results_ac2/" containing a plethora of analysis results. 
  
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
    
  - Execute all jobs until one or more specific rules: (until implies including)
    
    ```bash
    asscom2 --until roary abricate
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


![dag](https://github.com/cmkobel/assemblycomparator2/assets/5913696/3164f060-3b36-4d51-8cf7-29a50d87ec84)



**Hint:** Use `asscom2 --until <rulename> [<rulename2>...]` to run one or several specific analyses only. The rulename for each analysis to pick is listed below:

### For each assembly
  - `copy` [any2fasta](https://github.com/tseemann/any2fasta) Wide input format support and validation.
  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) Lengths and GC-content of individual contigs.
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) Generic assembly statistics.
  - `busco` [Busco](https://busco.ezlab.org/) Estimate assembly completeness and contamination.
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) Estimate assembly completeness and contamination.
  - `prokka` [prokka](https://github.com/tseemann/prokka) Genomic annotation.
  - `diamond_kegg` [diamond](https://github.com/bbuchfink/diamond) Run prokka-called proteins through the checkm2 database (uniref100 with KEGG-KOs).
  - `kegg_pathway` [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) KEGG ortholog-based pathway enrichment analysis.
  - `dbcan` [dbCAN4](https://github.com/linnabrown/run_dbcan) Annotation of carbohydrate-active "CAZyme" enzymes (lacking in report).
  - `interproscan` [InterProScan](https://github.com/ebi-pf-team/interproscan) Protein function using Tigrfam, Hamap and Pfam.
  - `abricate` [abricate](https://github.com/tseemann/abricate) Virulence and resistance gene identification.
  - `mlst` [mlst](https://github.com/tseemann/mlst) Multi locus sequence typing.
  - `kraken2` [kraken2](https://ccb.jhu.edu/software/kraken2/) Species identification.
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) Species recognition.
  

### For each group
  - `roary` [roary](https://sanger-pathogens.github.io/Roary/) (pan and core genome)
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) (core genome pairwise snp-distances)
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) (phylogenetic tree of the core genome)
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) (super fast distance measurement)
  - **A nice report easy to share with your friends ([demos](readme-demos.md))**


#### Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, kraken2, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, kraken2, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with `asscom2 --until meta` or `asscom2 --until isolate`



# Installation via Conda/Mamba

This is the recommended way to install Assemblycomparator2. Please note that it is highly recommended that you have apptainer on your system as it makes Assemblycomparator2 able to use a compressed docker-image that speeds up installation significantly. If you are not able to [install Apptainer](https://apptainer.org/docs/user/main/quick_start.html#installation-request), Assemblycomparator2 will default to use Mamba to install all snakemake workflow rule environments separately.

First, you need to install a Conda-based Python3 distribution.
The recommended choice is [Miniforge](https://github.com/conda-forge/miniforge#install) which not only provides the required Python and Conda commands, 
but also includes Mamba an extremely fast and robust replacement for the Conda package manager which is highly recommended.
The default conda solver is a bit slow and sometimes has issues with [selecting the latest package releases](https://github.com/conda/conda/issues/9905). 
Therefore, we recommend to in any case use Mamba.

In case you don't use Mambaforge you can always install Mamba into any other Conda-based Python distribution with

```bash

conda install -n base -c conda-forge mamba

```

Assemblycomparator2 can be installed into its own isolated environment by first creating the environment and then installing.

```bash

mamba create -c conda-forge -c bioconda -n asscom2 assemblycomparator2

```

Installing into isolated environments is best practice in order to avoid side effects with other packages.

 


## Optionally: Testing the installation

Now you will be able to run asscom2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run asscom2 it will show the message "Pulling singularity image docker://cmkobel/assemblycomparator2." This might take some time depending on your network bandwidth as it downloads a +4GB docker image that contains all the conda environments needed for each analysis.
```bash
conda activate asscom2

# First, create an empty directory and enter.
mkdir test_ac2_install
cd test_ac2_install

# Copy some test metagenomic assemblies from the test directory.
cp $CONDA_PREFIX/assemblycomparator2/tests/MAGs/*.fasta .

# Should take about a minute to complete the "fast" pseudo-rule.
asscom2 --until fast

# You can then investigate the report document that has been generated.
# open results_ac2/report_test_ac2_install.html

# Downloads all databases (~ 200 GB).
asscom2 --until downloads

# Run the full pipeline (~ 1 cpu-hour per genome).
asscom2
```




<table><tr><td>

#### Shared database

If you are working on a shared computational resource like a laboratory workstation or a HPC you might want to share a database directory so that each user will not have to redundantly download each database. To set this up, the first user must decide on a directory and set reading and writing permissions for the group of users that should be able to use the database. Writing permissions are necessary for the "database representative" flags that snakemake uses to keep track of the presence of the databases. Setting this custom path is a matter of defining the "ASSCOM2_DATABASES" environment variable. You can put this into your ~/.bashrc or execute the command before using asscom2.

```bash
export ASSCOM2_DATABASES="/absolute/path/to/shared_databases/asscom2_v2.5.8+"
```


</td></tr></table>



<table><tr><td>

#### Rule development

If you want to develop new [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-and-rules) in the assemblycomparator2 pipeline, you should consider following the [development version installation instructions](readme-development.md). The development version is purely conda-based so you can effect the next version apptainer-compatible docker image. 

</td></tr></table>







## Future functionality 

In the future we might add some of the following software packages into assemblycomparator2.

**Sample basis**
  - [Bakta](https://github.com/oschwengers/bakta/) (Rapid & standardized annotation of bacterial genomes, MAGs & plasmids)
  - MCF (Module Completion Factor as an extension to the KEGG-based pathway enrichment analysis that helps interpretation) 
  - [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) (Identify possible replication origins, and thereby help identify chromids)
  - [RFplasmid](https://github.com/aldertzomer/RFPlasmid) (Identify plasmids using the pentamer-random-forest method)
  - [Kaptive](https://github.com/katholt/Kaptive) (surface polysaccharide loci for Klebsiella and Acinetobacter baumannii) 
  - [mash screen](https://mash.readthedocs.io/en/latest/tutorials.html) (recognition of plasmids-of-interest)


**Batch basis**

  - [IQ-tree](http://www.iqtree.org/) (phylogenetic tree of core genome with bootstrapping)
  - GC3-profiling ("fingerprinting" of the distribution of GC-content)
  - Identification of horizontally transferred genes?
  - Recombination in core genome using the Bruen's PHI statistic
  - [GenAPI](https://github.com/MigleSur/GenAPI) (alternative to roary)



## Citation

If you use assemblycomparator2, you can support further funding by bumping the citation count on this one:

  - Kobel CM. *assemblycomparator2* **GitHub** https://github.com/cmkobel/assemblycomparator2

assemblycomparator2 would not have existed, if it hadn't been for the integrated software packages and their databases.

  - W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962. 
  - Richard Challis. (2017). rjchallis/assembly-stats 17.02 (17.02). Zenodo. https://doi.org/10.5281/zenodo.322347
  - Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
  - Chklovski, A., Parks, D.H., Woodcroft, B.J. et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods 20, 1203–1212 (2023). https://doi.org/10.1038/s41592-023-01940-w
  - Seemann T., Prokka: rapid prokaryotic genome annotation, Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063
  - Yin Y*, Mao X*, Yang JC, Chen X, Mao F and Xu Y, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Res. 2012 Jul;40(Web Server issue):W445-51 
  - P. Jones et al., Bioinformatics (2014), PMID: 24451626
  - Seemann T, Abricate, Github https://github.com/tseemann/abricate
  - Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018 Sep 24;3:124. doi: 10.12688/wellcomeopenres.14826.1. PMID: 30345391; PMCID: PMC6192448.
  - Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
  - Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics. 2019 Nov 15;36(6):1925–7. doi: 10.1093/bioinformatics/btz848. Epub ahead of print. PMID: 31730192; PMCID: PMC7703759.
  - Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

---


Development is active and will continue.

Assemblycomparator2 genomes to report pipeline. Copyright (C) 2019-2023 [Carl M. Kobel](https://github.com/cmkobel) [GNU GPL v3](https://github.com/cmkobel/assemblycomparator2/blob/master/LICENSE)
  
  
