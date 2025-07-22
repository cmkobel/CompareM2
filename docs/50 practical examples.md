# Practical examples

!!! note
    (This page is under construction. Will be finished before the end of august 2025.)

CompareM2 contains different tools to cover varying aspects of microbial genome analysis. Here we showcase some different ways to use CompareM2 to achieve insights into the functional potentials of microbial genomes, and how they relate to one another.


## Tutorial 1 - Publicly available *Enterococcus faecium*

The quickest way to get started with CompareM2 is to use publicly available genomes. Here, we'll investigate clusters of the *Enterococcus faecium* species, which are human commensals and pathogens, sometimes carrying vancomycin resistance.

First, we'll gain a quick phylogenetic overview over some publicly available genomes. A sufficient tree can be computed with mashtree.

After following the [installation instructions](https://comparem2.readthedocs.io/en/latest/10%20installation/), the following can be run on the command line

```bash
conda activate comparem2_launcher

mkdir E._faecium_public && cd E._faecium_public

comparem2 --config add_refseq=GCF...,GCF... --until mashtree
```

## Tutorial 2 - Metagenomic characterization

Assemblies of microbial genomes come in many different forms. A common way to investigate the diversity of ecological niches is to assemble and bin the microbes in a specific sample. By assembling a metagenome and binning it into discrete species representatives, it is possible to paint a picture of the functional potential they ?direct on their environment.

Here, we'll use genomes from the ?? study, which can be downloaded here.


{!resources/footer.md!}
