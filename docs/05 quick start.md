# Quick start ⚡️

If you have a set of bacterial or archaeal genomes, either MAGs or isolates, that you wish to analyze and compare—CompareM2 is the tool for you. 

## 1) Install

Assuming you already have [conda/mamba installed](https://github.com/conda-forge/miniforge#install), CompareM2 can be installed in a single step.

```bash

mamba create -c conda-forge -c bioconda -n comparem2 comparem2

```

## 2) Run

CompareM2 has a large number of available tools ([list](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#for-each-sample)), here we just want to run the ones that are fast.

```bash

comparem2 --config input_genomes="path/to/my/genomes_*.fna" --until fast 

```

## 3) Explore

When CompareM2 is done running, you can start exploring the results. The first thing to do, is to explore the dynamic report that shows the most important results from each analysis.

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


If you wish to run the [full](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/) rulegraph, simply remove the _until_-parameter when running CompareM2.

If you have any problems using CompareM2, you're very welcome to file an issue on the git repo: [https://github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues).


{!resources/footer.md!}

