

# What analyses does it do?

## The CompareM2 pipeline

CompareM2 can perform a large number of bioinformatic tasks which are each defined in "rules". Rules can be thought of as task templates that generate the computer code that performs the same computational task on many different input files (e.g. inputted genomes). Some of the rules are interdependent on one another and form chains of results that successively feed into more high level analyses down the line - hence pipeline. Below is a visualization of the dependency relationships between these rules.
 

![dag2 pdf](https://github.com/user-attachments/assets/855674a4-d80b-4892-8b14-5d87ad7de86b)

The visualized pipeline is a directed acyclic graph. The direction goes from start (copy) in the top where inputted genomes are copied into the results directory, to the end (all) in the bottom where results are collected before being directed into the dynamic report document. 


In many cases it may be desirable to run only a subset of the analyses on a given set of microbial assemblies. CompareM2 uses Snakemake's powerful pipeline execution to direct the order of execution for the many rules and their dependencies. This means that the user can select to run only specific parts of the rule graph. This can be achieved in practice by using the "until" parameter. 


By defining a rule to run *until*, the Snakemake executor can select to use only the parts of the pipeline, that lead up to the production of the desired analysis. One example is to only compute a core-pan genome using panaroo. To do this, first the inputted genome must be copied into the results directory (rule all), then the genome is annotated (rule annotate) using bakta og prokka, and then finally the core-pan genome can be computed across several samples (rule panaroo). To run this exact chain of dependencies of panaroo including itself, the user can simply run `comparem2 --until panaroo`. By running this command, only the relevant jobs for generating the core-pan genome are computed.


!!! note "TL;DR"
    Use `comparem2 --until <rule> [<another rule>...]` to run one or several specific analyses only. The rule names for each analysis to pick is listed below:

## Included analyses

Below is a comprehensive list of all rules (analyses) available in CompareM2.

### For each sample

First, independent analyses are run on each of the input genomic assembly files.

  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) Lengths and GC-content of individual contigs.
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) Generic assembly statistics.
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) Assembly completeness and contamination.
  - `bakta` [bakta](https://github.com/oschwengers/bakta) Genomic annotation of Bacteria (lacking in report, but used downstream by other tools).
  - `prokka` [prokka](https://github.com/tseemann/prokka) Legacy genomic annotation of Archaea and Bacteria. 
  - `kegg_pathway` [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) KEGG ortholog-based pathway enrichment analysis.
  - `dbcan` [dbCAN4](https://github.com/linnabrown/run_dbcan) Annotation of carbohydrate-active enzymes (CAZymes) (lacking in report).
  - `antismash` [antismash](https://docs.antismash.secondarymetabolites.org/) Detection of biosynthesis genes (lacking in report).
  - `eggnog` [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/) Functional annotation.
  - `interproscan` [InterProScan](https://github.com/ebi-pf-team/interproscan) Protein function using Tigrfam, Hamap and Pfam (lacking in report).
  - `amrfinder` [AMRFinderPlus](https://github.com/ncbi/amr/) Virulence and resistance gene identification.
  - `mlst` [mlst](https://github.com/tseemann/mlst) Multi locus sequence typing.
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) Species recognition.
  - `gapseq` [Gapseq](https://gapseq.readthedocs.io/en/latest/) Genome scale metabolic reconstruction gapfilling and modeling.
  

## Across samples

Then, on the basis of the analysis of each input genomic assembly, these analyses are run across all samples.

  - `panaroo` [panaroo](https://github.com/gtonkinhill/panaroo) Pan and core genome.
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) Pairwise snp-distances on the core genome.
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) Phylogenetic tree of the core genome.
  - `iqtree` [IQ-tree](http://www.iqtree.org/) Phylogenetic Tree of core genome with bootstrapping (lacking in report).
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) Super fast distance measurement and neighbor joining.
  - `treecluster` [TreeCluster](https://github.com/niemasd/TreeCluster) Clustering of phylogenetic trees (lacking in report).
  - **A nice report easy to share with your friends (See demos [below](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report))**


## Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, mlst, abricate, assembly-stats, gtdbtk, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, gtdbtk, checkm2, mashtree.


!!! note "Hint"
    You can run one of these pseudorules just like any other rulename with `comparem2 --until meta` or `comparem2 --until isolate`



## Rendered report

These demo reports are available for inspiration while you wait for your own report to complete.

  - [report_strachan_campylo.html](https://github.com/cmkobel/comparem2/raw/master/tests/strachan_campylo/report_strachan_campylo.html.zip)

    32 Campylobacter genomes, Metagenome and genome sequencing from the rumen epithelial wall of dairy cattle. From Nature 2022 - Strachan et al. (doi.<nolink />org/10.1038/s41564-022-01300-y).
    
  - [report_Methanoflorens.html](https://github.com/cmkobel/comparem2/raw/master/tests/Methanoflorens/report_Methanoflorens.html.zip)
  
    6 Methanoflorens (archaeal) genomes. Representants of Bog-38 which are part of GTDB.
    
  - [report_]
  
    500 Lysobacteraceae representants. 



{!resources/footer.md!}
