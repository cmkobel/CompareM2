# CompareM2

[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/cmkobel/comparem2/latest-fast.yaml)](https://github.com/cmkobel/CompareM2/actions) [![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/comparem2?label=Bioconda%20downloads&color=%2300CC00)](https://comparem2.readthedocs.io/en/latest/10%20installation/) [![Docker Pulls](https://img.shields.io/docker/pulls/cmkobel/comparem2?label=docker%20pulls)](https://comparem2.readthedocs.io/en/latest/10%20installation/) [![Documentation Status](https://readthedocs.org/projects/comparem2/badge/?version=latest)](https://comparem2.readthedocs.io/en/latest/?badge=latest) [![Conda Version](https://img.shields.io/conda/v/bioconda/comparem2)](https://anaconda.org/bioconda/comparem2) [![https://doi.org/10.1093/bioinformatics/btaf517](https://img.shields.io/badge/doi%20%28OUP%29-10.1093%2Fbioinformatics%2Fbtaf517-blue.svg)](https://doi.org/10.1093/bioinformatics/btaf517)

> **Note:** Looking for the original CompareM (AAI and codon usage)? See [github.com/donovan-h-parks/CompareM](https://github.com/donovan-h-parks/CompareM).

CompareM2 is a genomes-to-report pipeline for comparative analysis of bacterial and archaeal genomes. It takes genome assemblies — isolates or MAGs — and runs 30+ analysis tools, producing a portable HTML report with publication-ready graphics.

<a href="https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report"><img width="220" style="width: 220px" alt="report document logo" align="right" src="https://github.com/cmkobel/comparem2/assets/5913696/e5f9b72c-2137-4850-8779-a5528d8ccbaf"></a>

## Features

- **30+ integrated analyses** — quality control, annotation, functional annotation, phylogenetics, pan/core genomes, AMR profiling, metabolic modeling, and more. [Full list](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/).
- **Dynamic HTML report** — collects central results into a single portable file with interpretable figures and text. Adapts to partial runs — only includes sections for completed analyses. [See examples](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report).
- **Assembly-agnostic** — works strictly downstream of assembly and binning. Accepts genomes from any sequencing technology or source.
- **Scalable** — runs on local workstations (recommended >= 64 GiB RAM) or HPC clusters (SLURM/PBS). Scales approximately linearly with input size thanks to Snakemake's parallel job scheduling.
- **Easy to install** — single command via pixi or mamba. All dependencies are managed via conda environments or a pre-built Docker/Apptainer image.
- **Configurable** — [passthrough arguments](https://comparem2.readthedocs.io/en/latest/20%20usage/#passthrough-arguments) forward any parameter to any underlying tool. Add NCBI reference genomes by accession.

## Quick start

```bash
# Install
pixi global install -c conda-forge -c bioconda comparem2

# Run fast analyses on your genomes
comparem2 --config input_genomes="*.fna" --until fast

# Run the full pipeline
comparem2
```

See the [documentation](https://comparem2.readthedocs.io) for installation options, usage details, and the full list of analyses.

## Citation

Kobel C.M., Aho V.T.E., Øyås O., Nørskov-Lauritsen N., Woodcroft B.J., Pope P.B. CompareM2 is a genomes-to-report pipeline for comparing microbial genomes. *Bioinformatics* 41(9), btaf517 (2025). [doi:10.1093/bioinformatics/btaf517](https://doi.org/10.1093/bioinformatics/btaf517)

## Links

- **Documentation**: [comparem2.readthedocs.io](https://comparem2.readthedocs.io)
- **Issues and questions**: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues)
