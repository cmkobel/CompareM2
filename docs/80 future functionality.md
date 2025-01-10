

# Future functionality 

In the future we might add some of the following software packages into CompareM2. This document serves as a backlog of tools that we want to integrate when time allows.

**Assembly basis (within each sample)**
  
  - [AlphaFold](https://github.com/google-deepmind/alphafold) Neural network protein folding prediction genome annotation.
  - Integration of the [DRAM](https://github.com/WrightonLabCSU/DRAM) databases for easier metabolic interpretation.
  - [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html) Identification of possible replication origins of chromids.
  - [RFplasmid](https://github.com/aldertzomer/RFPlasmid) Identification of plasmids using the pentamer-random-forest method.
  - [gapseq](https://github.com/jotech/gapseq/tree/master) GEMs, pathway completeness and much more.
  - [distillR](https://github.com/anttonalberdi/distillR) High level functional annotation using graph based metabolic capacity indices.


**Batch basis (across all samples)**

  - GC3-profiling "fingerprinting" of the distribution of GC-content.
  - Recombination in core genome using the Bruen's PHI statistic or ClonalFrameML. (requires stable synteny calculation)
  - Identification of horizontally transferred genes?

Please [add an issue on the repository](https://github.com/cmkobel/comparem2/issues) if you have any ideas or requests.


{!resources/footer.md!}
