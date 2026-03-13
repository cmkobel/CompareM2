

# What analyses does it do?

## Pipeline overview

CompareM2 runs a directed acyclic graph (DAG) of interdependent analysis rules. Each rule performs a specific bioinformatic task on one or more input genomes. Rules form chains where outputs feed into downstream analyses.

![dag2 pdf](https://github.com/user-attachments/assets/855674a4-d80b-4892-8b14-5d87ad7de86b)

[\[open figure\]](https://github.com/user-attachments/assets/855674a4-d80b-4892-8b14-5d87ad7de86b)

The pipeline flows from top (copying input genomes) to bottom (collecting results into the report).

### N-dependent output selection

CompareM2 automatically selects which analyses to run based on the number of input genomes (N):

| N | Analyses enabled |
|---|---|
| 0 | Database downloads only |
| 1+ | Per-genome analyses: annotation, QC, functional annotation, metabolic modeling |
| 2+ | Pairwise comparisons: panaroo, snp-dists, mashtree, treecluster |
| 3+ | Phylogenetics: fasttree, iqtree, bootstrap mashtree |

### Running specific analyses

Use `--until` to run only specific rules and their dependencies:

```bash
comparem2 --until panaroo
```

This runs only what is needed to produce the panaroo output: copy, annotate, and panaroo itself.

!!! note
    Use `comparem2 --until <rule> [<rule2>...]` to run one or more specific analyses. Rule names are listed below.


## Included analyses

### Quality control

| Rule | Tool | Description |
|---|---|---|
| `assembly_stats` | [Assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | Assembly statistics (N50, total length, etc.) |
| `sequence_lengths` | [SeqKit](https://bioinf.shenwei.me/seqkit/usage/) | Per-contig lengths and GC content |
| `checkm2` | [CheckM2](https://github.com/chklovski/CheckM2/) | Completeness and contamination estimation |

### Annotation

The annotation output is used by many downstream tools.

| Rule | Tool | Description |
|---|---|---|
| `bakta` | [Bakta](https://github.com/oschwengers/bakta) | Genomic annotation of bacteria (default annotator) |
| `prokka` | [Prokka](https://github.com/tseemann/prokka) | Genomic annotation of bacteria and archaea |

### Functional annotation

| Rule | Tool | Description |
|---|---|---|
| `eggnog` | [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper/) | Functional annotation via orthology |
| `interproscan` | [InterProScan](https://github.com/ebi-pf-team/interproscan) | Protein function (TIGRFAM, Hamap, Pfam) |
| `dbcan` | [dbCAN](https://github.com/linnabrown/run_dbcan) | Carbohydrate-active enzyme (CAZyme) annotation |
| `kegg_pathway` | [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/) | KEGG pathway enrichment analysis |
| `amrfinder` | [AMRFinderPlus](https://github.com/ncbi/amr/) | Antimicrobial resistance and virulence genes |
| `mlst` | [MLST](https://github.com/tseemann/mlst) | Multi-locus sequence typing |
| `gapseq_find` / `gapseq_fill` | [gapseq](https://gapseq.readthedocs.io/en/latest/) | Metabolic pathway prediction and gap-filling |
| `antismash` | [antiSMASH](https://docs.antismash.secondarymetabolites.org/) | Biosynthetic gene cluster detection |

### Core and pan genomes

| Rule | Tool | Description |
|---|---|---|
| `panaroo` | [Panaroo](https://github.com/gtonkinhill/panaroo) | Pan and core genome analysis |

### Phylogenetics and taxonomy

| Rule | Tool | Description |
|---|---|---|
| `mashtree` | [Mashtree](https://github.com/lskatz/mashtree) | Fast distance estimation and neighbor-joining tree |
| `bootstrap_mashtree` | [Mashtree](https://github.com/lskatz/mashtree) | Mashtree with bootstrap support (N >= 3) |
| `fasttree` | [FastTree](http://www.microbesonline.org/fasttree/) | Core genome phylogeny |
| `iqtree` | [IQ-TREE](http://www.iqtree.org/) | Core genome phylogeny with bootstrapping |
| `gtdbtk` | [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/) | Taxonomic classification |
| `snp_dists` | [snp-dists](https://github.com/tseemann/snp-dists) | Pairwise SNP distances on the core genome |
| `treecluster` | [TreeCluster](https://github.com/niemasd/TreeCluster) | Phylogenetic tree clustering |

### Dynamic report

The report is always generated and collects results from all completed analyses.

  - `report` — A portable HTML report with interpretable results and publication-ready graphics. See [demo reports below](#rendered-report).


## Pseudo-rules

Pseudo-rules are shortcuts to run curated subsets of the pipeline:

| Pseudo-rule | Included analyses |
|---|---|
| `fast` | sequence_lengths, assembly-stats, mashtree |
| `meta` | annotation, assembly-stats, sequence_lengths, checkm2, eggnog, kegg_pathway, dbcan, interproscan, gtdbtk, mashtree |
| `isolate` | annotation, assembly-stats, sequence_lengths, eggnog, kegg_pathway, gtdbtk, mlst, amrfinder, panaroo, fasttree, snp-dists, mashtree |
| `downloads` | All database download rules |
| `report` | Re-render the report |

!!! note "Hint"
    Run a pseudo-rule like any other rule: `comparem2 --until meta` or `comparem2 --until isolate`


## Rendered report

These demo reports are available for download:

  - [report_strachan_campylo.html](https://github.com/cmkobel/comparem2/raw/master/tests/strachan_campylo/report_strachan_campylo.html.zip) — 32 *Campylobacter* genomes from Strachan et al. (Nature 2022, doi.org/10.1038/s41564-022-01300-y). Metagenome and genome sequencing from the rumen epithelial wall of dairy cattle.

  - [report_Methanoflorens.html](https://github.com/cmkobel/comparem2/raw/master/tests/Methanoflorens/report_Methanoflorens.html.zip) — 6 *Methanoflorens* (archaeal) genomes. Representatives of Bog-38 from GTDB.



{!resources/footer.md!}
