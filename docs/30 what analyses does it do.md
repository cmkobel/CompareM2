

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

#### `assembly_stats` — [Assembly-stats](https://github.com/sanger-pathogens/assembly-stats)

Computes assembly statistics including N50 (the length of the smallest contig that, together with the longer contigs, covers at least half of the genome), total assembly length, number of contigs, and other summary metrics.

**Requires:** N ≥ 1. No database download needed.

---

#### `sequence_lengths` — [SeqKit](https://bioinf.shenwei.me/seqkit/usage/)

Extracts per-contig lengths and GC content from each input genome. The report visualizes each contig as a bar colored by GC content, giving a quick overview of assembly fragmentation and composition bias.

**Requires:** N ≥ 1. No database download needed.

---

#### `checkm2` — [CheckM2](https://github.com/chklovski/CheckM2/)

Estimates genome completeness and contamination using machine learning on a universal gene set. Essential for assessing the quality of metagenome-assembled genomes (MAGs).

**Requires:** N ≥ 1. Downloads the CheckM2 DIAMOND database (~3.5 GB) on first run.

---

### Annotation

The annotation output is used by many downstream tools (eggNOG, dbCAN, InterProScan, Panaroo, etc.). Choose one annotator via the `annotator` config setting. NCBI-sourced genomes automatically use their bundled annotation instead.

#### `bakta` — [Bakta](https://github.com/oschwengers/bakta) (default)

Rapid and standardized annotation of bacterial genomes. Bakta uses a comprehensive pre-built database and produces consistent locus tags suitable for comparative analyses.

**Requires:** N ≥ 1. Downloads the Bakta database on first run.

| Parameter | Default | Description |
|---|---|---|
| `set_bakta--translation-table` | `11` | Genetic code translation table |
| `set_bakta--gram` | `"?"` | Gram type for signal peptide prediction (`+`, `-`, or `?`) |
| `set_bakta--meta` | *(unset)* | Enable metagenome mode (flag) |

---

#### `prokka` — [Prokka](https://github.com/tseemann/prokka)

Whole genome annotation for bacteria and archaea. Alternative to Bakta — select with `annotator: "prokka"` in config.

**Requires:** N ≥ 1. No database download needed.

| Parameter | Default | Description |
|---|---|---|
| `set_prokka--compliant` | *(flag set)* | Force Genbank/ENA/DDJB compliance |
| `set_prokka--kingdom` | `bacteria` | Annotation kingdom (`archaea`, `bacteria`, `mitochondria`, `viruses`) |
| `set_prokka--gram` | *(unset)* | Gram type (`neg`, `pos`) |
| `set_prokka--rfam` | *(unset)* | Enable Rfam search for ncRNAs (flag) |

---

### Functional annotation

#### `eggnog` — [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper/)

Functional annotation through orthology assignment. Maps predicted proteins against a database of pre-computed orthologous groups to transfer functional information including COG categories, KEGG orthologs, and Gene Ontology terms.

**Requires:** N ≥ 1. Downloads the eggNOG database on first run.

| Parameter | Default | Description |
|---|---|---|
| `set_eggnog-m` | `diamond` | Search mode (`diamond`, `mmseqs`, `hmmer`) |

---

#### `interproscan` — [InterProScan](https://github.com/ebi-pf-team/interproscan)

Classifies proteins into families and predicts domains and important sites using multiple member databases (TIGRFAM, Hamap, Pfam by default).

**Requires:** N ≥ 1. No database download needed (InterProScan bundles its own data).

| Parameter | Default | Description |
|---|---|---|
| `set_interproscan--applications` | `TIGRFAM,Hamap,Pfam` | Comma-separated list of member databases to run |
| `set_interproscan--goterms` | *(flag set)* | Include Gene Ontology terms |
| `set_interproscan--pathways` | *(flag set)* | Include pathway annotations |

---

#### `dbcan` — [dbCAN](https://github.com/linnabrown/run_dbcan)

Annotates carbohydrate-active enzymes (CAZymes) by searching against the dbCAN HMM and DIAMOND databases. Useful for studying carbohydrate metabolism, degradation, and biosynthesis capabilities.

**Requires:** N ≥ 1. Downloads the dbCAN database on first run.

---

#### `kegg_pathway` — [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/)

KEGG pathway enrichment analysis. Predicted proteins are searched against the UniRef100-KO database (≥85% coverage, ≥50% identity), and [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)'s `enricher` function computes Benjamini-Hochberg adjusted p-values for pathway enrichment per genome.

**Requires:** N ≥ 1. Uses the CheckM2 DIAMOND database (shared download).

---

#### `amrfinder` — [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)

Identifies antimicrobial resistance genes, point mutations, and virulence and stress resistance genes in assembled nucleotide and protein sequences using NCBI's curated reference database.

**Requires:** N ≥ 1. No separate database download needed.

---

#### `mlst` — [MLST](https://github.com/tseemann/mlst)

Multi-locus sequence typing using the [PubMLST](https://pubmlst.org/) database. Automatically detects the appropriate MLST scheme for each genome and assigns a sequence type.

**Requires:** N ≥ 1. No separate database download needed.

| Parameter | Default | Description |
|---|---|---|
| `set_mlst--scheme` | *(auto-detect)* | Force a specific MLST scheme (e.g., `efaecium`, `saureus`) |

---

#### `gapseq_find` / `gapseq_fill` — [gapseq](https://gapseq.readthedocs.io/en/latest/)

Predicts metabolic pathways (`gapseq_find`) and reconstructs gap-filled genome-scale metabolic models (`gapseq_fill`). The two-step process first identifies pathways and transporters, then fills gaps in the metabolic network to produce a functional model.

**Requires:** N ≥ 1. No separate database download needed.

| Parameter | Default | Description |
|---|---|---|
| `set_gapseq_find-t` | `auto` | Taxonomic range for reference sequences (`Bacteria`, `Archaea`, `auto`) |
| `set_gapseq_fill_draft-b` | `auto` | Biomass reaction to use |

---

#### `antismash` — [antiSMASH](https://docs.antismash.secondarymetabolites.org/)

Detects and characterizes biosynthetic gene clusters (BGCs) for secondary metabolites including antibiotics, siderophores, and terpenes.

**Requires:** N ≥ 1. Downloads the antiSMASH database on first run.

---

#### `carveme` — [CarveMe](https://github.com/cdanielmachado/carveme)

Automated reconstruction of genome-scale metabolic models from annotated genomes. Produces SBML models suitable for flux balance analysis.

**Requires:** N ≥ 1. No separate database download needed.

| Parameter | Default | Description |
|---|---|---|
| `set_carveme--gapfill` | `LB` | Growth media for gap-filling (`M9`, `LB`, or comma-separated) |
| `set_carveme--solver` | `scip` | LP solver to use |

---

### Core and pan genomes

#### `panaroo` — [Panaroo](https://gtonkinhill.github.io/panaroo/#/)

Computes the pan and core genome across input genomes. The core genome contains genes conserved across all (or nearly all) samples, while the pan genome is the union of all genes. Panaroo also produces a core genome alignment used by downstream phylogenetic tools.

**Requires:** N ≥ 2.

| Parameter | Default | Description |
|---|---|---|
| `set_panaroo--clean-mode` | `sensitive` | Error-correction stringency (`strict`, `moderate`, `sensitive`) |
| `set_panaroo--core_threshold` | `0.95` | Fraction of samples a gene must appear in to be considered "core" |
| `set_panaroo--threshold` | `0.98` | Sequence identity threshold for clustering |
| `set_panaroo-a` | `core` | Alignment output type (`core`, `pan`) |
| `set_panaroo-f` | `0.7` | Protein family sequence identity threshold |
| `set_panaroo--remove-invalid-genes` | *(flag set)* | Exclude genes with unusual length or premature stop codons |

---

### Phylogenetics and taxonomy

#### `mashtree` — [Mashtree](https://github.com/lskatz/mashtree)

Computes an approximation of ANI using the MinHash distance measure and builds a neighbor-joining tree. Fast enough for hundreds of genomes. The resulting tree is unrooted.

**Requires:** N ≥ 2.

| Parameter | Default | Description |
|---|---|---|
| `set_mashtree--genomesize` | `5000000` | Expected genome size (bp) |
| `set_mashtree--mindepth` | `5` | Minimum k-mer depth |
| `set_mashtree--kmerlength` | `21` | K-mer length |
| `set_mashtree--sketch-size` | `10000` | Sketch size for MinHash |

---

#### `bootstrap_mashtree` — [Mashtree](https://github.com/lskatz/mashtree)

Mashtree with bootstrap support values. Inherits all `set_mashtree` parameters from above.

**Requires:** N ≥ 3.

| Parameter | Default | Description |
|---|---|---|
| `set_bootstrap_mashtree--reps` | `100` | Number of bootstrap replicates |

---

#### `fasttree` — [FastTree](http://www.microbesonline.org/fasttree/)

Builds an approximately-maximum-likelihood phylogenetic tree from the core genome alignment produced by Panaroo. Faster than IQ-TREE but with less rigorous statistical support.

**Requires:** N ≥ 3. Depends on Panaroo core genome alignment.

| Parameter | Default | Description |
|---|---|---|
| `set_fasttree-gtr` | *(flag set)* | Use the generalized time-reversible (GTR) model |

---

#### `iqtree` — [IQ-TREE](http://www.iqtree.org/)

Maximum-likelihood phylogenetic inference with bootstrap support from the core genome alignment. More thorough than FastTree, providing formal model selection and statistical branch support.

**Requires:** N ≥ 3. Depends on Panaroo core genome alignment.

| Parameter | Default | Description |
|---|---|---|
| `set_iqtree--boot` | `100` | Number of bootstrap replicates |
| `set_iqtree-m` | `GTR` | Substitution model |

---

#### `gtdbtk` — [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/)

Taxonomic classification using the Genome Taxonomy Database (GTDB). Assigns species names by measuring average nucleotide identity (ANI) and relative evolutionary divergence (RED) against reference sequences.

**Requires:** N ≥ 1. Downloads the GTDB database (~85 GB) on first run.

| Parameter | Default | Description |
|---|---|---|
| `set_gtdbtk--keep_intermediates` | *(flag set)* | Retain intermediate files |

---

#### `snp_dists` — [snp-dists](https://github.com/tseemann/snp-dists)

Counts pairwise SNP differences on the core genome alignment. Note: SNP distances are not adjusted for transition/transversion bias and give a ballpark indication of divergence rather than a true evolutionary distance. Highly sensitive to the core/pan genome size ratio.

**Requires:** N ≥ 2. Depends on Panaroo core genome alignment.

---

#### `treecluster` — [TreeCluster](https://github.com/niemasd/TreeCluster)

Clusters genomes on a phylogenetic tree using a distance threshold. Useful for defining operational taxonomic units or outbreak clusters.

**Requires:** N ≥ 2. Runs on the Mashtree output.

| Parameter | Default | Description |
|---|---|---|
| `set_treecluster--method` | `max_clade` | Clustering method (see [TreeCluster docs](https://github.com/niemasd/TreeCluster) for all options) |
| `set_treecluster--threshold` | `0.05` | Distance threshold for cluster assignment |

---

### Dynamic report

The report is always generated and collects results from all completed analyses. Only sections for tools that ran successfully are included.

  - `report` — A portable HTML report with interpretable results and publication-ready graphics. See [demo reports below](#rendered-report).


## Passthrough parameters

Any tool parameter can be forwarded from the CompareM2 config using the `set_` prefix. The naming convention is `set_<tool><flag>: <value>`. Flag-only arguments (no value) use an empty string `""`.

For example, to change the IQ-TREE substitution model and increase bootstrap replicates:

```yaml
# In config/config.yaml or via --config
set_iqtree-m: GTR+G
set_iqtree--boot: 1000
```

Or on the command line:

```bash
comparem2 --config set_iqtree-m=GTR+G set_iqtree--boot=1000
```

To add a flag argument (no value), set it to an empty string:

```bash
comparem2 --config 'set_prokka--rfam=""'
```

To remove a default passthrough parameter, delete (or comment out) the corresponding line in `config/config.yaml`. Default values for all passthrough parameters are listed in the tool sections above.

!!! note
    Check each tool's own documentation (linked above) for the full list of available flags.


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
