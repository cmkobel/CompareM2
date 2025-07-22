# Tutorials

CompareM2 contains many very different tools to cover many aspects of microbial genome analysis. Here we showcase some different ways to use CompareM2 to achieve insights into the functional potentials of microbial genomes, and how they relate to one another.


## Tutorial 1 - Publicly available *Enterococcus faecium*

The quickest way to get started with CompareM2 is to use publicly available genomes. Here, we'll investigate clusters of the *Enterococcus faecium* species, which are human commensals and pathogens, sometimes carrying vancomycin resistance.

First, we'll gain a quick phylogenetic overview over some publicly available genomes. A sufficient tree can be computed with mashtree.

```
mkdir E._faecium_public && cd E._faecium_public

```bash
comparem2 --config add_refseq=GCF...,GCF... --until mashtree
```


{!resources/footer.md!}
