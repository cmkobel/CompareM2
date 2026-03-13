# CompareM2

CompareM2 is a genomes-to-report pipeline for comparing microbial genomes. Given bacterial or archaeal genome assemblies — isolates or MAGs — it runs 30+ analysis tools and produces a portable HTML report with publication-ready graphics.

!!! note
    If you're looking for the original CompareM (AAI and codon usage), see [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM).

## What it does

- **Quality control** — assembly statistics, completeness and contamination (CheckM2)
- **Annotation** — Bakta or Prokka, with functional annotation via eggNOG, InterProScan, dbCAN, and antiSMASH
- **Resistance and virulence** — AMRFinderPlus, MLST
- **Phylogenetics and taxonomy** — Mashtree, FastTree, IQ-TREE, GTDB-Tk, SNP distances
- **Pan/core genomes** — Panaroo
- **Metabolic modeling** — gapseq, KEGG pathway enrichment

See the [full list of analyses](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/).

## Get started

```bash
# Install
pixi global install -c conda-forge -c bioconda comparem2

# Run fast analyses
comparem2 --config input_genomes="*.fna" --until fast

# Run everything
comparem2
```

See the [quick start guide](https://comparem2.readthedocs.io/en/latest/05%20quick%20start/) or [installation instructions](https://comparem2.readthedocs.io/en/latest/10%20installation/) for details.

## How it works

CompareM2 is a Snakemake pipeline. It automatically selects which analyses to run based on the number of input genomes, manages all software dependencies via conda environments or a pre-built Docker/Apptainer image, and collects results into a single HTML report.

It is **assembly-agnostic** — it works strictly downstream of assembly and binning. Bring genomes from any sequencing technology or source, and CompareM2 handles the rest. This is a deliberate design choice: read mapping, assembly, and binning are highly dependent on sequencing technology and are best handled by specialized tools. CompareM2 focuses on what comes after.

The dynamic report only includes sections for analyses that completed, so it adapts to partial runs. It is designed to be interpretable by non-bioinformaticians, with explanatory text and figures alongside the results.

Benchmarking showed that CompareM2 scales approximately linearly with input size thanks to Snakemake's parallel job scheduling, and is significantly faster than comparable tools like Tormes and Bactopia ([Kobel et al. 2025](https://doi.org/10.1093/bioinformatics/btaf517)).

## Links

- **Documentation**: [comparem2.readthedocs.io](https://comparem2.readthedocs.io)
- **Source code**: [github.com/cmkobel/CompareM2](https://github.com/cmkobel/CompareM2)
- **Issues and questions**: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues)
- **Citation**: Kobel et al. *Bioinformatics* 41, btaf517 (2025). [doi:10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

{!resources/footer.md!}
