# Quick start

## Install

Assuming you already have conda/mamba installed.

```bash

mamba create -c conda-forge -c bioconda -n comparem2 comparem2

```

## Run

Using `--until fast` to only run the fast analyses.

```bash

comparem2 --config input_genomes="path/to/my/genomes_*.fna" --until fast 

open results_comparem2/report_*.html

```

If you wish to run the [full](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/) rulegraph, simply remove the _until_-parameter.




{!resources/footer.md!}

