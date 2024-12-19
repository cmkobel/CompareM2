# Quick start ⚡️

You have a set of bacterial or archaeal genomes, either MAGs or isolates, that you wish to analyze and compare.

## 1) Install

Assuming you already have conda/mamba installed.

```bash

mamba create -c conda-forge -c bioconda -n comparem2 comparem2

```

## 2) Run

Using `--until fast` to only run the fast analyses.

```bash

comparem2 --config input_genomes="path/to/my/genomes_*.fna" --until fast 

```

## 3) Explore

Open results_comparem2/report_*.html ([demo](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report)) in your browser to gain a quick oversight of the results.

The full results can be explored from the results_comparem2/ directory:
 
```bash

tree results_comparem2/ -L 1
#> results_comparem2/
#> ├── amrfinder/
#> ├── assembly-stats/
#> ├── benchmarks/
#> ├── checkm2/
#> ├── fasttree/
#> ├── gtdbtk/
#> ├── iqtree/
#> ├── kegg_pathway/
#> ├── mashtree/
#> ├── metadata.tsv
#> ├── mlst/
#> ├── panaroo/
#> ├── report_<title>.html
#> ├── samples/
#> │  └── <sample>/
#> │     ├── antismash/
#> │     ├── bakta/
#> │     ├── dbcan/
#> │     ├── eggnog/
#> │     ├── <sample>.fna
#> │     ├── interproscan/
#> │     ├── prokka/
#> │     └── sequence_lengths/
#> ├── snp-dists/
#> ├── tables/
#> ├── treecluster/
#> └── version_info.txt

```


If you wish to run the [full](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/) rulegraph, simply remove the _until_-parameter.




{!resources/footer.md!}
