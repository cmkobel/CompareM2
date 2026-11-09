# Future functionality

## What CompareM2 deliberately does not do

Some things are out of scope by choice and not for lack of time, and knowing
which is which saves everyone an issue.

  - **Deep functional and metabolic interpretation of a single genome.** This
    is where [DRAM](https://github.com/WrightonLabCSU/DRAM) and
    [nf-core/funcscan](https://nf-co.re/funcscan) are strong, and they
    complement CompareM2 rather than overlapping with it.
  - **Read-level processing.** Mapping, assembly and binning are highly
    dependent on sequencing technology. CompareM2 works strictly downstream of
    them, which is what makes it technology-independent.

Left out of v3, each for a reason of its own: eggNOG-mapper, InterProScan and
dbCAN, which are the deep functional annotation above; gapseq, because CarveMe
builds the same kind of model in minutes rather than hours; Prokka, because
Bakta is the one annotator and weighs about half as much; clusterProfiler, an R
package in a pipeline that no longer has R; and IQ-TREE, with FastTree building
the core-genome tree.

## Candidates

Genuinely wanted, not yet in:

  - **A faster provisional taxonomy** than GTDB-Tk, whose 60.8 GB database is
    97% of the measured database total. Less pressing since r232 more than
    halved it. The requirement is that it takes *assemblies*. sylph was tried
    and removed because it profiles metagenomic reads, which is not the
    question being asked here.
  - **Recombination detection** in the core genome, with Bruen's PHI statistic
    or ClonalFrameML. The core-genome tree currently assumes one shared history
    for all genes, which recombination breaks.
  - **Plasmid identification**, with
    [RFplasmid](https://github.com/aldertzomer/RFPlasmid) or similar. Panaroo's
    `--clean-mode strict` can remove rare plasmids, so this is partly a
    correction.
  - **Replication origin identification**, with
    [Oriloc](https://doi.org/10.1093/bioinformatics/16.6.560), which ships in
    the R package `seqinr`.
  - **GC3 profiling**, synonymous GC-content fingerprinting.
  - **Horizontally transferred gene identification.** The pangenome section can
    already point at candidates: a gene pattern shared by genomes that are not
    neighbours in the tree.
  - **Biosynthetic gene clusters**, with antiSMASH. It was selected for v3 and
    then dropped, because it pins `biopython 1.78` and `diamond 2.1.11` against
    newer versions other tools here need, and at the time every tool had to
    share one environment. They no longer do. A dependency set that conflicts
    with the others is the kind of reason a seventh environment needs, so the
    conflict is no longer what keeps it out. Its databases would have to be
    measured first.

Anything added declares its database's download size, measured from
`content-length` where the database has one and declared unmeasured where it
does not, so that the cost printed before a run stays honest. It joins one of
the six conda environments, unless it brings a dependency ecosystem of its own,
which is a reason for a seventh, written down with the measurement that
justified it. An environment per tool is not. See `DESIGN.md`.

Please [open an issue](https://github.com/cmkobel/CompareM2/issues) if you have
ideas or requests.
