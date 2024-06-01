

# What analyses does it do?

Below is the graph the shows the order of execution of all possible analyses "rules" in Assemblycomparator2:


![dag](https://github.com/cmkobel/assemblycomparator2/assets/5913696/db6b58ac-cec0-43fe-b06f-18048ef3b642)


This figure does not show the pseudo rules such as `meta`, `isolate`, `fast`, etc.

**Hint:** Use `asscom2 --until <rule> [<another rule>...]` to run one or several specific analyses only. The rule names for each analysis to pick is listed below:

## For each sample

First, independent analyses are run on each of the input genomic assembly files.

  - `copy` [any2fasta](https://github.com/tseemann/any2fasta) Wide input format support and validation.
  - `sequence_lengths` [seqkit](https://bioinf.shenwei.me/seqkit/usage/) Lengths and GC-content of individual contigs.
  - `assembly_stats` [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) Generic assembly statistics.
  - `busco` [Busco](https://busco.ezlab.org/) Estimate assembly completeness and contamination.
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) Estimate assembly completeness and contamination.
  - `prokka` [prokka](https://github.com/tseemann/prokka) Genomic annotation of Archaea and Bacteria. 
  - `bakta` [bakta](https://github.com/oschwengers/bakta) Genomic annotation of Bacteria (lacking in report, but used downstream by other tools).
  - `kegg_pathway` [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) KEGG ortholog-based pathway enrichment analysis.
  - `dbcan` [dbCAN4](https://github.com/linnabrown/run_dbcan) Annotation of carbohydrate-active "CAZyme" enzymes (lacking in report).
  - `antismash` [antismash](https://docs.antismash.secondarymetabolites.org/) Detection of biosynthesis genes (lacking in report).
  - `eggnog` [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/) Functional annotation.
  - `interproscan` [InterProScan](https://github.com/ebi-pf-team/interproscan) Protein function using Tigrfam, Hamap and Pfam (lacking in report).
  - `abricate` [abricate](https://github.com/tseemann/abricate) Virulence and resistance gene identification.
  - `mlst` [mlst](https://github.com/tseemann/mlst) Multi locus sequence typing.
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) Species recognition.
  

## Across samples

Then on the basis of the analysis of each input genomic assembly, these analyses are run across all samples.

  - `panaroo` [panaroo](https://github.com/gtonkinhill/panaroo) Pan and core genome.
  - `snp_dists` [snp-dists](https://github.com/tseemann/snp-dists) Core genome pairwise snp-distances.
  - `fasttree` [FastTree](http://www.microbesonline.org/fasttree/) Phylogenetic tree of the core genome.
  - `iqtree` [IQ-tree](http://www.iqtree.org/) Phylogenetic Tree of core genome with bootstrapping (lacking in report).
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) Super fast distance measurement
  - `treecluster` [TreeCluster](https://github.com/niemasd/TreeCluster) Clustering of phylogenetic trees (lacking in report).
  - **A nice report easy to share with your friends (See demos [below](https://assemblycomparator2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#demo-reports))**


## Pseudo-rules

There are also a few pseudo targets defined. For instance `fast` which runs sequence_lengths, assembly-stats and mashtree. There is also one named `isolate` which runs all the analyses that are relevant for clinical isolates (sequence_lengths, prokka, mlst, abricate, assembly-stats, gtdbtk, busco, checkm2, roary, snp-dists, fasttree, mashtree) as well as one named `meta` which runs the analyses that are relevant to metagenomes (aka. MAGs), these are sequence_lengths, prokka, gtdbtk, busco, checkm2, mashtree.

**Hint:** You can run one of these pseudorules just like any other rulename with `asscom2 --until meta` or `asscom2 --until isolate`



## Rendered report

These demo reports are available for inspiration while you wait for your own report to complete.

  - [report_strachan_campylo.html](https://github.com/cmkobel/assemblycomparator2/raw/master/tests/strachan_campylo/report_strachan_campylo.html.zip)

    32 Campylobacter genomes, Metagenome and genome sequencing from the rumen epithelial wall of dairy cattle. From Nature 2022 - Strachan et al. (doi.<nolink />org/10.1038/s41564-022-01300-y).
  - 





