

# What analyses does it do?

## The CompareM2 pipeline

CompareM2 can perform a large number of bioinformatic tasks which are each defined in "rules". Rules can be thought of as task templates that generate the computer code that performs the same computational task on many different input files (e.g. inputted genomes). Some of the rules are interdependent on one another and form chains of results that successively feed into more high level analyses down the line - hence pipeline. Below is a visualization of the dependency relationships between these rules. The arrows show the order of the individual rules, when the pipeline runs.
 

![dag2 pdf](https://github.com/user-attachments/assets/855674a4-d80b-4892-8b14-5d87ad7de86b)

[\[open figure\]](https://github.com/user-attachments/assets/855674a4-d80b-4892-8b14-5d87ad7de86b)

The visualized pipeline is a directed acyclic graph. The direction goes from start (copy) in the top where inputted genomes are copied into the results directory, to the end (all) in the bottom where results are collected before being directed into the dynamic report document. 


In many cases it may be desirable to run only a subset of the analyses on a given set of microbial assemblies. CompareM2 uses Snakemake's powerful pipeline execution to direct the order of execution for the many rules and their dependencies. This means that the user can select to run only specific parts of the rule graph. This can be achieved in practice by using the "until" parameter. 


### Controlling execution with the `--until` parameter

By defining a rule to run *until*, the Snakemake executor can select to use only the parts of the pipeline, that lead up to the production of the desired analysis. One example is to only compute a core-pan genome using panaroo. To do this, first the inputted genome must be copied into the results directory (rule all), then the genome is annotated (rule annotate) using bakta og prokka, and then finally the core-pan genome can be computed across several samples (rule panaroo). To run this exact chain of dependencies of panaroo including itself, the user can simply run `comparem2 --until panaroo`. By running this command, only the relevant jobs for generating the core-pan genome are computed.


!!! note Summary
    Use `comparem2 --until <rule> [<another rule>...]` to run one or several specific analyses only. The rule names for each analysis to pick is listed in the next section:

## Included analyses

Below is a comprehensive list of all rules (analyses) available in CompareM2.
                         
    
### Quality Control
  - `assembly_stats` [Assembly-stats](https://github.com/sanger-pathogens/assembly-stats) generic assembly statistics.
  - `sequence_lengths` [Seqkit](https://bioinf.shenwei.me/seqkit/usage/) lengths and GC-content of individual contigs.
  - `checkm2` [CheckM2](https://github.com/chklovski/CheckM2/) assembly completeness and contamination.

### Annotation

The output of the annotation tools is used downstream in the other tools. 

  - `bakta` [Bakta](https://github.com/oschwengers/bakta) genomic annotation of Bacteria (lacking in report, but used downstream by other tools).
  - `prokka` [Prokka](https://github.com/tseemann/prokka) legacy genomic annotation of Archaea and Bacteria. 
  
### Advanced annotation
  - `eggnog` [Eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/) functional annotation.
  - `interproscan` [Interproscan](https://github.com/ebi-pf-team/interproscan) protein function using Tigrfam, Hamap and Pfam (lacking in report).
  - `dbcan` [Dbcan](https://github.com/linnabrown/run_dbcan) annotation of carbohydrate-active enzymes (CAZymes) (lacking in report).
  - `kegg_pathway` [Clusterprofiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) KEGG ortholog-based pathway enrichment analysis.
  - `amrfinder` [Amrfinderplus](https://github.com/ncbi/amr/) virulence and resistance gene identification.
  - `mlst` [Mlst](https://github.com/tseemann/mlst) multi locus sequence typing.
  - `gapseq` [Gapseq](https://gapseq.readthedocs.io/en/latest/) genome scale metabolic reconstruction gapfilling and modeling.
  - `antismash` [Antismash](https://docs.antismash.secondarymetabolites.org/) detection of biosynthesis genes (lacking in report).
  
### Core-pan genomes
  - `panaroo` [Panaroo](https://github.com/gtonkinhill/panaroo) pan and core genomes.

### Phylogenetic or Taxonomic
  - `mashtree` [Mashtree](https://github.com/lskatz/mashtree) super fast distance measurement and neighbor joining.
  - `fasttree` [Fasttree](http://www.microbesonline.org/fasttree/) phylogenetic tree of the core genome.
  - `gtdbtk` [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) species identification.
  - `iqtree` [Iq-tree](http://www.iqtree.org/) phylogenetic tree of core genome with bootstrapping (lacking in report).
  - `treecluster` [Treecluster](https://github.com/niemasd/TreeCluster) clustering of phylogenetic trees (lacking in report).
  - `snp_dists` [Snp-dists](https://github.com/tseemann/snp-dists) pairwise snp-distances on the core genome.
  
### Dynamic report

The report always runs and collects all the results.
  
  - `report` **A nice report easy to share with your friends (See demos [below](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report))**
  


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
