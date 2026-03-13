# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CompareM2 is a Snakemake-based bioinformatics pipeline for comparing microbial genomes (bacteria and archaea, including MAGs). It is **not** a Python package — it is a workflow application with a Python launcher script. It only runs on Linux. It runs 36+ analysis tools and produces a publication-ready HTML report.

## Architecture

Three-layer design:
1. **Launcher** (`./comparem2`) — Python script that sets environment variables (`COMPAREM2_BASE`, `COMPAREM2_PROFILE`, `COMPAREM2_DATABASES`), auto-detects Apptainer vs Conda, parses CLI args, and invokes Snakemake
2. **Snakemake pipeline** (`workflow/Snakefile` + `workflow/rules/*.smk`) — 11 rule files organized as `sample_*.smk` (per-genome) and `batch_*.smk` (cross-genome), plus `downloads.smk` and `gather.smk`
3. **Dynamic report** (`dynamic_report/`) — separate sub-pipeline generating the final HTML report

Key config: `config/config.yaml`. Passthrough parameters use `set_` prefix to forward tool-specific arguments (e.g., `set_iqtree--boot: 100`).

Conda environments for each tool live in `workflow/envs/*.yaml` (25 YAML files). Execution profiles (conda vs apptainer, local vs HPC) are in `profile/`.

**Execution modes:** Local conda, HPC (SLURM/PBS profiles in `profile/`), or containerized (Docker/Apptainer).

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

## Testing & CI/CD

No unit test framework — testing is done via CI/CD (`.github/workflows/`) running the pipeline on test datasets in `tests/`. Key workflows:
- `latest-dry-run.yaml` — dry-run validation
- `latest-fast.yaml` — weekly fast test (conda)
- `latest-full.yaml` — full test suite
- `stable-conda.yaml`, `stable-apptainer_A/B.yaml` — stability/release tests

Test datasets: `tests/E._faecium/`, `tests/E._faecium_plasmids/`, `tests/MAGs/`, `tests/Methanoflorens/`, `tests/strachan_campylo/`, `tests/nocore/`

## Key Design Patterns

- **Annotator selection**: `bakta` (default) or `prokka`, configured via `annotator` in config. NCBI-sourced genomes automatically use their bundled annotation. Downstream rules consume whichever annotation via symlinked `.annotation/` directory.
- **NCBI integration**: Reference genomes can be added via `add_ncbi` config parameter (comma-separated accessions). Downloads are cached in `.ncbi_cache/`.
- **Passthrough parameters**: Any tool parameter can be forwarded using `set_<tool><flag>: <value>` in config. Flag-only arguments use empty string.
- **Origin tracking**: Local genomes vs NCBI accessions follow different processing paths but merge for downstream analysis.
- **Metadata**: Managed via pandas DataFrames in the Snakefile for sample tracking and status reporting.

## Configuration

Main config: `config/config.yaml`. Key settings:
- `input_genomes` — glob pattern for input files
- `fofn` — file-of-filenames (alternative to glob)
- `annotator` — "bakta" (default) or "prokka"
- `output_directory` — default "results_comparem2"
- ~28 default `set_<tool>` passthrough parameters (users can add more)
- `title` — custom report title (defaults to working directory name)

**Pseudo targets** (use with `--until`): `fast`, `meta`, `isolate`, `downloads`, `report`
