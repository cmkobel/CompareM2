# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CompareM2 is a Snakemake-based bioinformatics pipeline for comparing microbial genomes (bacteria and archaea, including MAGs). It runs 36+ analysis tools and produces a publication-ready HTML report.

## Architecture

Three-layer design:
1. **Launcher** (`./comparem2`) — Python script that sets environment variables (`COMPAREM2_BASE`, `COMPAREM2_PROFILE`, `COMPAREM2_DATABASES`), parses CLI args, and invokes snakemake
2. **Snakemake pipeline** (`workflow/Snakefile` + `workflow/rules/*.smk`) — 11 rule files organized as `sample_*.smk` (per-genome) and `batch_*.smk` (cross-genome), plus `downloads.smk` and `gather.smk`
3. **Dynamic report** (`dynamic_report/`) — separate sub-pipeline generating the final HTML report

Key config: `config/config.yaml`. Passthrough parameters use `set_` prefix to forward tool-specific arguments (e.g., `set_iqtree--boot: 100`).

Conda environments for each tool live in `workflow/envs/*.yaml`. Execution profiles (conda vs apptainer, local vs HPC) are in `profile/`.

## Common Commands

```bash
# Set up development environment
mamba env create -y -f environment.yaml -n comparem2_dev
conda activate comparem2_dev

# Development testing (conda)
export COMPAREM2_BASE="$(realpath ~/comparem2)"
export COMPAREM2_PROFILE="${COMPAREM2_BASE}/profile/conda/default"
${COMPAREM2_BASE}/comparem2 --config input_genomes="${COMPAREM2_BASE}/tests/E._faecium/*.fna" --until fast

# Same but with apptainer
export COMPAREM2_PROFILE="${COMPAREM2_BASE}/profile/apptainer/default"

# Dry run (validate pipeline structure)
./comparem2 --config input_genomes="*.fna" --dry-run

# Run up to a specific rule
./comparem2 --until <rule_name>

# Show pipeline status
./comparem2 --status

# Download databases
./comparem2 --downloads

# Generate DAG visualization
conda install anaconda::graphviz
./comparem2 --forceall --rulegraph | dot -Tpng > dag.png

# Update Dockerfile
touch dummy1.fa dummy2.fa dummy3.fa
./comparem2 --config add_ncbi=GCF2987 --containerize | grep -A 10000 "FROM condaforge" | grep -B 10000 "mamba clean --all -y" > Dockerfile
```

## Versioning

Version must be bumped in three places simultaneously:
1. `comparem2` (line ~7, `__version__`)
2. `workflow/Snakefile` (version variable)
3. `changelog.txt`

Docker image is pinned to minor version (e.g., `docker://cmkobel/comparem2:v2.16`). Database directory is also namespaced by minor version (`cm2_v2.16/`).

## CI/CD

GitHub Actions workflows in `.github/workflows/`:
- `latest-dry-run.yaml` — structure validation
- `latest-fast.yaml` — quick functional test
- `latest-full.yaml` — complete pipeline test
- `stable-conda.yaml`, `stable-apptainer_A/B.yaml` — release testing

Test datasets in `tests/` (E. faecium, MAGs, Methanoflorens, Campylobacter).

## Key Design Patterns

- **Annotator selection**: `bakta` (default) or `prokka`, configured via `annotator` in config. Downstream rules consume whichever annotator output is selected.
- **NCBI integration**: Reference genomes can be added via `add_ncbi` config parameter; handled through checkpoint-based conditional logic in the Snakefile.
- **Passthrough parameters**: Any tool parameter can be forwarded using `set_<tool><flag>: <value>` in config. Flag-only arguments use empty string.
- **Metadata**: Managed via pandas DataFrames in the Snakefile for sample tracking and status reporting.
