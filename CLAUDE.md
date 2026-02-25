# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CompareM2 is a Snakemake-based bioinformatics pipeline for microbial genome analysis. It is **not** a Python package — it is a workflow application with a Python launcher script. It only runs on Linux.

## Architecture

**Two-stage pipeline:**
1. **Main workflow** (`workflow/Snakefile` + rules in `workflow/rules/`): Snakemake DAG executing parallel genome analyses
2. **Dynamic report** (`dynamic_report/`): Separate Snakemake subpipeline generating HTML reports from R Markdown sections

**Key entry point:** `./comparem2` (Python launcher, 373 lines) — sets environment variables (`COMPAREM2_BASE`, `COMPAREM2_PROFILE`, `COMPAREM2_DATABASES`), auto-detects Apptainer vs Conda, then invokes Snakemake.

**Rule organization** (`workflow/rules/`):
- `sample_*.smk` — per-genome analyses (annotation, QC, clinical/AMR)
- `batch_*.smk` — cross-genome analyses (phylogeny, pan/core genome, batch QC)
- `downloads.smk` — database management
- `gather.smk` — output aggregation

**Execution modes:** Local conda, HPC (SLURM/PBS profiles in `profile/`), or containerized (Docker/Apptainer). Profiles live in `profile/conda/` and `profile/apptainer/`.

**Key design patterns:**
- **Passthrough parameters:** `set_<rule>--parameter=value` in config allows per-tool customization without code changes
- **Annotator abstraction:** Bakta (default), Prokka, or NCBI annotation with transparent switching via symlinks
- **Origin tracking:** Local genomes vs NCBI accessions follow different processing paths but merge for downstream analysis

## Versioning

Version must be bumped in three places (kept in sync):
1. `comparem2` line 7 (`__version__`)
2. `workflow/Snakefile` line 7 (`__version__`)
3. `changelog.txt`

Docker images are tagged by minor version (e.g., `cmkobel/comparem2:v2.16`).

## Development Setup

```bash
cd ~
git clone https://github.com/cmkobel/comparem2.git comparem2
cd comparem2
mamba env create -y -f environment.yaml -n comparem2_dev
conda activate comparem2_dev
# Force conda over Apptainer (if Apptainer is installed):
export COMPAREM2_PROFILE="$(realpath profile/conda/default)"
```

## Common Commands

```bash
# Dry run (validate pipeline without executing)
./comparem2 --config input_genomes="*.fna" --until fast --dry-run

# Run fast rule set on test data
unzip tests/E._faecium/fna.zip
./comparem2 --config input_genomes="*.fna" --until fast

# Check pipeline status
./comparem2 --status

# Download databases
./comparem2 --downloads

# Update Dockerfile from Snakemake containerize output
touch dummy1.fa dummy2.fa dummy3.fa
./comparem2 --config add_ncbi=GCF2987 --containerize | grep -A 10000 "FROM condaforge" | grep -B 10000 "mamba clean --all -y" > Dockerfile

# Update DAG diagram
conda install anaconda::graphviz
comparem2 --forceall --rulegraph | dot -Tpng > dag.png
```

## Testing

No unit test framework — testing is done via CI/CD (`.github/workflows/`) running the pipeline on test datasets in `tests/`. Key workflows:
- `latest-fast.yaml` — weekly fast test (conda)
- `latest-full.yaml` — full test suite
- `latest-dry-run.yaml` — dry-run validation
- `stable-conda.yaml`, `stable-apptainer_A/B.yaml` — stability tests

Test datasets: `tests/E._faecium/`, `tests/MAGs/`, `tests/Methanoflorens/`, `tests/strachan_campylo/`

## Configuration

Main config: `config/config.yaml`. Key settings:
- `input_genomes` — glob pattern for input files
- `fofn` — file-of-filenames (alternative to glob)
- `annotator` — "bakta" (default) or "prokka"
- `output_directory` — default "results_comparem2"
- 50+ `set_<tool>` passthrough parameters

Tool-specific conda environments: `workflow/envs/` (28 YAML files).
