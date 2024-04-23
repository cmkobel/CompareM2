# Assemblycomparator2

[![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/assemblycomparator2?label=docker%20pulls)](#installation)  [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/assemblycomparator2?label=Bioconda%20downloads&color=%2300CC00)](#installation) [![conda build](https://img.shields.io/conda/v/bioconda/assemblycomparator2)](#installation)

🧬 Assemblycomparator2 "asscom2" is a genomes-to-report pipeline. It accepts prokaryotic genomic assemblies and compares them in many different ways. 

🦠 Being designed to analyze assemblies of both isolates and metagenomes (MAGs), it is useful for anyone working with microbial genomics.

💾 [Installing](#installation) Assemblycomparator2 on your system gives you access to [15](https://github.com/cmkobel/assemblycomparator2#what-analyses-does-it-do) powerful state-of-the-art tools for analysis of prokaryotic genomes which will accelerate your research. It is easy to use and can be used by non-bioinformaticians.

<img alt="asscom2 animation" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/623f6b42-2de6-457c-8f0d-3b3e5d646967">



👩‍🔬 Assemblycomparator2 integrates several analyses that yield scientific results about genomic assemblies on several levels: Quality control, annotation, function and species calling as well as comparative analyses like computation of core/pan genomes and phylogenetics. 

<img width="150" alt="snakemake logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/7188e748-9d37-43ae-a5d5-100e9560df1f">

🐍 Assemblycomparator2 works by calling a Snakemake workflow that can be easily modified to use different parameters for the  underlying tools.

<a href="https://github.com/cmkobel/assemblycomparator2/blob/master/readme-demos.md"><img height="192" alt="report document logo" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/e5f9b72c-2137-4850-8779-a5528d8ccbaf"></a>

📙 Central results are dynamically integrated in a compact portable report .html-document. It can be browsed in any web browser and can be easily shared as a single file. This report is generated even if some jobs in the pipeline fail. See [examples](readme-demos.md).

🧑‍💻 Assemblycomparator2 can be run either on a local workstation (recommended >= 64GiB RAM), or a HPC (high performance computing) cluster. Both  Apptainer/Singularity/Docker images and conda environment definitions are available for all dependent software to run.

Assemblycomparator2 will -depending on circumstances- be renamed to "Genomisk" or something completely different.

## Usage examples

Make a directory with the prokaryotic genomic assembly-files -or metagenomic bins- you want to investigate with Assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
Assemblycomparator2 will then create a sub-directory, named "results_ac2/" containing a plethora of analysis results. 
  
  - Execute a "dry run". That is, to show what will be run without actually doing it.

    ```bash
    asscom2 --dry-run
    ```

  - Run Assemblycomparator2 on the genomes in the current directory:

    ```bash
    asscom2
    ```
    

##### A bit more advanced controls 

  - Run analyses that are relevant to metagenomic assemblies only (as opposed to isolates):

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

Below is the graph the shows the order of execution of all possible analyses "rules" in Assemblycomparator2:


![dag](https://github.com/cmkobel/assemblycomparator2/assets/5913696/3164f060-3b36-4d51-8cf7-29a50d87ec84)

This figure does not show the pseudo rules such as `meta`, `isolate`, `fast`, etc.

**Hint:** Use `asscom2 --until <rulename> [<rulename2>...]` to run one or several specific analyses only. The rulename for each analysis to pick is listed below:

### For each sample

First, independent analyses are run on each of the input genomic assembly files.

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
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) Species recognition.
  

### Across samples

Then on the basis of the analysis of each input genomic assembly, these analyses are run across all samples.

  - `roary` [roary](https://sanger-pathogens.github.io/Roary/) Pan and core genome.
  - `motulizer` and `motupan` [mOTUlizer](https://github.com/moritzbuck/mOTUlizer) Analyze core-pan spectrum genome and gene clusters (lacking in report).
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) Core genome pairwise snp-distances.
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) Phylogenetic tree of the core genome.
  - `iqtree` [IQ-tree](http://www.iqtree.org/) Phylogenetic Tree of core genome with bootstrapping (lacking in report).
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) Super fast distance measurement
  - **A nice report easy to share with your friends ([demos](readme-demos.md))**


#### Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with `asscom2 --until meta` or `asscom2 --until isolate`



# Installation

<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/5b06b511-75c4-48cb-8ab8-f29b212ef6df">

It is highly recommended that you have [Apptainer](https://Apptainer.org/docs/user/main/quick_start.html#installation-request) on your system as it makes Assemblycomparator2 able to use a compressed Docker-image that speeds up installation significantly.



<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/c9d15678-b4a7-42be-b0de-b649479f6d74">

First, you need to install a Conda or Mamba package manager.
The recommended choice is [Miniforge](https://github.com/conda-forge/miniforge#install) which not only provides the required Python and Conda commands, 
but also includes Mamba - an extremely fast and robust replacement for the Conda package manager which is highly recommended.
<table><tr><td>
  
If you don't have Apptainer, Assemblycomparator2 will default to use Conda/Mamba to install all Snakemake workflow rule environments separately. In this case, you should set the conda channel priority to "flexible" with `conda config --set channel_priority flexible`

</td></tr></table>

<table><tr><td>
In case you don't use Miniforge you can always install Mamba into any other Conda-based Python distribution with:

```bash

conda install -n base -c conda-forge mamba

```

</td></tr></table>

<img width="150" align="right" src="https://github.com/cmkobel/assemblycomparator2/assets/5913696/6bc39697-7e90-49a0-a44e-64820f2c1024">

Finally, Assemblycomparator2 can be installed into its own environment with the correct channels like so:

```bash

mamba create -c conda-forge -c bioconda -n asscom2 assemblycomparator2

```

Installing into isolated environments is best practice in order to avoid side effects with other packages.


## Optionally: Testing the installation

Now you will be able to run asscom2. You can use the example data in path "tests/MAGs" to check that everything works. The first time you run asscom2 it will show the message "Pulling singularity image docker://cmkobel/assemblycomparator2." This might take some time depending on your network bandwidth as it downloads a +4GB Docker image that contains all the conda environments needed for each analysis.
```bash

# Activate the newly created conda environment containing the asscom2 launcher.
mamba activate asscom2

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



### Advanced configuration

#### Shared database

If you are working on a shared computational resource like a laboratory workstation or a HPC you might want to share a database directory so that each user will not have to redundantly download each database. To set this up, the first user must decide on a directory and set reading and writing permissions for the group of users that should be able to use the database. Writing permissions are necessary for the "database representative" flags that snakemake uses to keep track of the presence of the databases. Setting this custom path is a matter of defining the "ASSCOM2_DATABASES" environment variable. You can put this into your ~/.bashrc or execute the command before using asscom2.

```bash
export ASSCOM2_DATABASES="/absolute/path/to/shared_databases/asscom2_v2.5.8+"
```

#### HPC profiles for Snakemake

 If you have experience with snakemake and are working on a high performance computing cluster (HPC), you can modify and use the cluster configuration profiles in the "profiles/" directory. You can define the use of one of these profiles by setting the "ASSCOM2_PROFILE" environment variable. You can put this into your ~/.bashrc or execute the command before using asscom2. You can read more about snakemake profiles [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) or browse more default profiles [here](https://github.com/snakemake-profiles).

```bash
export ASSCOM2_PROFILE=${ASSCOM2_BASE}/profiles/apptainer/slurm-sigma2-saga
```
#### Rule development

If you want to develop new [rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-and-rules) in the Assemblycomparator2 pipeline, you should consider following the [development version installation instructions](readme-development.md). The development version is purely conda-based so you can affect the next version of the Apptainer-compatible Docker image. 







## Future functionality 

In the future we might add some of the following software packages into Assemblycomparator2.

**Assembly basis (within each sample)**
  - [Bakta](https://github.com/oschwengers/bakta/) Rapid & standardized annotation of bacterial genomes, MAGs & plasmids.
  - MCF Module Completion Fraction as an extension to the KEGG-based pathway enrichment analysis that helps interpretation.
  - [AlphaFold](https://github.com/google-deepmind/alphafold) Neural network protein folding prediction genome annotation.
  - Integration of the [DRAM](https://github.com/WrightonLabCSU/DRAM) databases for easier metabolic interpretation.
  - [Eggnogg-mapper](https://github.com/eggnogdb/eggnog-mapper) Functional annotation of novel sequences.
  - [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) Identification of possible replication origins of chromids.
  - [RFplasmid](https://github.com/aldertzomer/RFPlasmid) Identification of plasmids using the pentamer-random-forest method.
  - [mash screen](https://mash.readthedocs.io/en/latest/tutorials.html) Recognition of plasmids-of-interest.
  - [Kaptive](https://github.com/katholt/Kaptive) Identification of surface polysaccharide loci for Klebsiella and Acinetobacter baumannii.
  - [AMRFinderPlus](https://github.com/ncbi/amr/) Identification of AMR genes and their point mutations.


**Batch basis (across all samples)**

  - GC3-profiling "fingerprinting" of the distribution of GC-content.
  - Recombination in core genome using the Bruen's PHI statistic.
  - Identification of horizontally transferred genes?
  - [GenAPI](https://github.com/MigleSur/GenAPI) Roary alternative.



## Citation

If you use Assemblycomparator2, you can support further funding by bumping the citation count on this one:

  - Kobel CM. *Assemblycomparator2* **GitHub** https://github.com/cmkobel/assemblycomparator2

Assemblycomparator2 would not have existed, if it hadn't been for the integrated software packages and their databases. And please reach out to carl.mathias.kobel near nmbu.no if you think something is missing.

  - Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann S, Nahnsen S, Köster J. Sustainable data analysis with Snakemake. F1000Res. 2021 Apr 19;10:33. doi: 10.12688/f1000research.29032.2. PMCID: PMC8114187.
  - Seemann T, Goncalves da Silva A, Bulach DM, Schultz MB, Kwong JC, Howden BP. Nullarbor Github https://github.com/tseemann/nullarbor
  - W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962. 
  - Richard Challis. (2017). rjchallis/assembly-stats 17.02 (17.02). Zenodo. https://doi.org/10.5281/zenodo.322347
  - Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
  - Chklovski, A., Parks, D.H., Woodcroft, B.J. et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nat Methods 20, 1203–1212 (2023). https://doi.org/10.1038/s41592-023-01940-w
  - Seemann T., Prokka: rapid prokaryotic genome annotation, Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063
  - Yin Y*, Mao X*, Yang JC, Chen X, Mao F and Xu Y, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Res. 2012 Jul;40(Web Server issue):W445-51 
  - P. Jones et al., Bioinformatics (2014), PMID: 24451626
  - Seemann T, Abricate, Github https://github.com/tseemann/abricate
  - Katz et al., (2019). Mashtree: a rapid comparison of whole genome sequence files. Journal of Open Source Software, 4(44), 1762, https://doi.org/10.21105/joss.01762
  - Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018 Sep 24;3:124. doi: 10.12688/wellcomeopenres.14826.1. PMID: 30345391; PMCID: PMC6192448.
  - Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics. 2019 Nov 15;36(6):1925–7. doi: 10.1093/bioinformatics/btz848. Epub ahead of print. PMID: 31730192; PMCID: PMC7703759.
  - Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
  - B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534. https://doi.org/10.1093/molbev/msaa015
  - mOTUpan: a robust Bayesian approach to leverage metagenome assembled genomes for core-genome estimation Moritz Buck, Maliheh Mehrshad, and Stefan Bertilsson bioRxiv 2021.06.25.449606; doi: https://doi.org/10.1101/2021.06.25.449606

---


Development is active and will continue. We're actively looking for collaborators to join and synergize this project.

Assemblycomparator2 genomes to report pipeline. Copyright (C) 2024 [Carl M. Kobel](https://github.com/cmkobel) [GNU GPL v3](https://github.com/cmkobel/assemblycomparator2/blob/master/LICENSE)
  
  
