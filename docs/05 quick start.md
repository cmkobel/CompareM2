# Quick start

## Install

```bash

mamba create -c conda-forge -c bioconda -n comparem2 comparem2

```

## Run

Using `--until fast` to only run the fast analyses.

```bash

comparem2 --config input_genomes="path/to/my/genomes_*.fna" --until fast 

open results_comparem2/report_*.html

```




