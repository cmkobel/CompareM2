# Future functionality

## What CompareM2 deliberately does not do

CompareM2 is a wide view over many assemblies, not a deep view into each one.
Some
things are out of scope by choice rather than for lack of time, and knowing
which is which saves everyone an issue:

  - **Deep functional and metabolic interpretation** — this is where
    [DRAM2](https://github.com/WrightonLabCSU/DRAM) and
    [nf-core/funcscan](https://nf-co.re/funcscan) are already strong. CompareM2
    will not compete with them.
  - **Read-level processing** — mapping, assembly and binning are highly
    dependent on sequencing technology. CompareM2 works strictly downstream of
    them, which is what makes it technology-independent.
  - **Biosynthetic gene clusters** — antiSMASH was selected and then
    dropped: it pins `biopython 1.78` and `diamond 2.1.11` against the newer
    versions CheckM2, Bakta and GTDB-Tk need, and one solved environment was
    judged worth more than one BGC caller. funcscan ships four.

Dropped from v2 for the same reason: eggNOG-mapper, InterProScan, dbCAN,
gapseq, antiSMASH, IQ-TREE, clusterProfiler and Prokka.

## Candidates

Genuinely wanted, not yet in:

  - **A faster provisional taxonomy** than GTDB-Tk, whose 60.8 GB database is
    97% of the measured database total. Less pressing since r232 more than
    halved it. The requirement is that it takes *assemblies* —
    sylph was tried and removed because it profiles metagenomic reads, which is
    not the question being asked here.
  - **Aligned fraction alongside ANI.** skani computes it, but `triangle
    --full-matrix` emits identity only, and an ANI is reported once alignment
    covers as little as ~15% of a genome. Identity without coverage is half the
    picture.
  - **Recombination detection** in the core genome — Bruen's PHI statistic or
    ClonalFrameML. The core-genome tree currently assumes one shared history for
    all genes, which recombination breaks.
  - **Plasmid identification** — [RFplasmid](https://github.com/aldertzomer/RFPlasmid)
    or similar. Panaroo's `--clean-mode strict` can remove rare plasmids, so
    this is partly a correction.
  - **Replication origin identification** —
    [Oriloc](http://pbil.univ-lyon1.fr/software/Oriloc/oriloc.html).
  - **GC3 profiling** — synonymous GC-content fingerprinting.
  - **Horizontally transferred gene identification.** The pangenome section can
    already point at candidates: a gene pattern shared by genomes that are not
    neighbours in the tree.

Anything added has to declare its database size in bytes and solve into the
existing environment. Both are hard constraints, not preferences — see
`DESIGN.md`.

Please [open an issue](https://github.com/cmkobel/CompareM2/issues) if you have
ideas or requests.
