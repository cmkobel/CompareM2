# Quick start

If you have bacterial or archaeal genomes — isolates or MAGs — that you want to analyze and compare, CompareM2 is the tool for you.

## 1) Install

With [pixi](https://pixi.sh) on a Linux machine:

```bash
pixi global install -c conda-forge -c bioconda comparem2
```

## 2) Run

CompareM2 includes many analysis tools ([full list](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/)). To run only the fast ones:

```bash
comparem2 --config input_genomes="path/to/my/genomes_*.fna" --until fast
```

## 3) Explore

When CompareM2 finishes, open the report in your browser:

```
results_comparem2/report_<title>.html
```

This dynamic HTML report summarizes the most important results from each analysis. See a [demo report](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#rendered-report).

The full results are organized in `results_comparem2/`:

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
#> ├── visuals/
#> ├── treecluster/
#> └── version_info.txt
```

To run the [full pipeline](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/), simply omit the `--until` parameter.

If you have any problems, please file an issue: [github.com/cmkobel/CompareM2/issues](https://github.com/cmkobel/CompareM2/issues).


{!resources/footer.md!}
