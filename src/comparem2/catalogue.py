"""The v3 tool set.

Fourteen tools, selected 2026-09-01. See DESIGN.md for what was dropped and why.

**Pin minimum versions.** An unpinned bioconda solve will happily pick a
years-old build to satisfy some other constraint, and the result resolves
cleanly but crashes at runtime. Both instances found so far were of exactly
this shape:

- `bakta` unpinned resolved to 1.8.1 (2023), which calls `pyrodigal.OrfFinder`
  — renamed `GeneFinder` in pyrodigal 3.x.
- `panaroo` unpinned resolved to 1.1.2 (2020), which imports `Bio.Alphabet`
  — removed in Biopython 1.78.

Neither is a solver failure. Only running the tool reveals them, which is why
the verification table in DESIGN.md tracks execution rather than installation.

Database sizes marked `size=` were measured from `content-length`; the rest are
`None` and must not be guessed at.
"""

from __future__ import annotations

from .tools import Context, Database, Registry, Scope, Tool

# --- Databases -----------------------------------------------------
# Sizes are measured, not estimated. `None` means nobody has measured it yet.

CHECKM2_DB = Database(
    name="checkm2",
    url="https://zenodo.org/records/14897628/files/checkm2_database.tar.gz?download=1",
    size=1_735_095_710,  # measured 2026-09-01
)

GTDB_DB = Database(
    name="gtdb",
    url=(
        "https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/"
        "auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz"
    ),
    size=141_442_235_198,  # measured 2026-09-01 — 91% of the whole install
)

# v2 used `--type full` (30 GB compressed / 84 GB on disk). v3 uses light.
# 1.3 GB is Bakta's documented figure, not a content-length measurement,
# because there is no static URL — `bakta_db download` picks the record.
#
# The database version is tied to the Bakta version: 1.12.x requires db 6.x and
# refuses 5.1 outright. So the download must run from the same environment as
# the tool, and pinning Bakta without re-fetching the database is a runtime
# failure rather than a solve failure.
BAKTA_DB = Database(name="bakta-light", url="bakta_db download --type light", size=None)

# Fetched by `amrfinder -u`; no static URL to measure. Small.
AMRFINDER_DB = Database(name="amrfinder", url="amrfinder -u", size=None)

# sylph needs a prebuilt GTDB sketch. Which release, and how large, is open.
SYLPH_DB = Database(name="sylph-gtdb", url="TBD", size=None)


# --- Quality and taxonomy ------------------------------------------

seqkit = Tool(
    name="seqkit",
    summary="Per-contig lengths and GC, and the assembly statistics derived from them.",
    scope=Scope.GENOME,
    conda=("bioconda::seqkit",),
    command=lambda c: [
        "seqkit", "fx2tab", "--name", "--length", "--gc",
        str(c.assembly), "-o", str(c.out("seqkit", "contigs.tsv")),
    ],
    outputs=lambda c: [c.out("seqkit", "contigs.tsv")],
)

checkm2 = Tool(
    name="checkm2",
    summary="Completeness and contamination for every genome.",
    scope=Scope.SET,
    conda=("bioconda::checkm2",),
    # --database_path takes the DIAMOND file itself, not the directory it sits
    # in; passing a directory fails with IsADirectoryError deep inside CheckM2.
    # The download step normalises the release's own name to checkm2.dmnd so
    # this path does not move when upstream renames the file.
    command=lambda c: [
        "checkm2", "predict", "--threads", str(c.threads),
        "--database_path", str(c.databases / "checkm2" / "checkm2.dmnd"),
        "--output-directory", str(c.out("checkm2")),
        "--input", *[str(a) for a in c.assemblies],
    ],
    outputs=lambda c: [c.out("checkm2", "quality_report.tsv")],
    database=CHECKM2_DB,
    threads=8,
    # CheckM2 pins an old DIAMOND (2.1.x) and cannot co-solve with a current
    # Bakta, which needs 2.2.x. Verified 2026-09-01: bakta>=1.10 solves with
    # all thirteen other tools, and fails only when checkm2 is added. Leaving
    # both unpinned "solves" by silently selecting bakta 1.8.1, which then
    # crashes on pyrodigal 3.x (OrfFinder was renamed GeneFinder). So this is
    # the one tool that gets its own environment.
    isolated=True,
)

gtdbtk = Tool(
    name="gtdbtk",
    summary="Taxonomic assignment against the GTDB reference tree.",
    scope=Scope.SET,
    conda=("bioconda::gtdbtk",),
    # Each genome lives in its own directory once inputs are canonicalised, so
    # --genome_dir cannot be used; the rule writes a batchfile first.
    # GTDB-Tk also emits bac120 and ar53 summaries separately, and the rule
    # concatenates them so archaea and bacteria appear in one table.
    command=lambda c: [
        "gtdbtk", "classify_wf", "--cpus", str(c.threads), "--skip_ani_screen",
        "--batchfile", str(c.out("gtdbtk", "batchfile.tsv")),
        "--out_dir", str(c.out("gtdbtk")),
    ],
    outputs=lambda c: [c.out("gtdbtk", "gtdbtk.summary.tsv")],
    database=GTDB_DB,
    threads=16,
)

sylph = Tool(
    name="sylph",
    summary="Fast provisional taxonomy, available long before GTDB-Tk finishes.",
    scope=Scope.SET,
    conda=("bioconda::sylph",),
    command=lambda c: [
        "sylph", "profile", "-t", str(c.threads),
        str(c.databases / "sylph" / "gtdb.syldb"),
        *[str(a) for a in c.assemblies],
        "-o", str(c.out("sylph", "profile.tsv")),
    ],
    outputs=lambda c: [c.out("sylph", "profile.tsv")],
    database=SYLPH_DB,
    threads=8,
)


# --- Annotation and screening --------------------------------------

bakta = Tool(
    name="bakta",
    summary="Structural and functional genome annotation.",
    scope=Scope.GENOME,
    conda=("bioconda::bakta>=1.10",),
    command=lambda c: [
        "bakta", "--db", str(c.databases / "bakta"), "--threads", str(c.threads),
        "--output", str(c.out("bakta")), "--prefix", str(c.sample), *c.args(),
        str(c.assembly),
    ],
    # Rules create their output directory before running, and bakta refuses to
    # write into an existing one. v2 passed --force for the same reason.
    params=(("--force", ""),),
    outputs=lambda c: [
        c.out("bakta", f"{c.sample}.gff3"),
        c.out("bakta", f"{c.sample}.faa"),
    ],
    database=BAKTA_DB,
    threads=8,
)

amrfinder = Tool(
    name="amrfinder",
    summary="Antimicrobial resistance and virulence genes.",
    scope=Scope.GENOME,
    conda=("bioconda::ncbi-amrfinderplus",),
    # Protein-only mode, as v2 used. Combined nucleotide+protein mode calls
    # better in principle, but AMRFinder cross-checks contig identifiers
    # between the GFF and the FASTA, and Bakta renames contigs to `contig_1`
    # while the input assembly keeps its original accessions — so combined mode
    # fails with `GFF contig id "contig_1" is not in the DNA FASTA file`.
    # The cost is losing point-mutation detection, which needs --organism.
    needs=("bakta",),
    command=lambda c: [
        "amrfinder", "-p", str(c.out("bakta", f"{c.sample}.faa")),
        "--plus", "--threads", str(c.threads), *c.args(),
    ],
    outputs=lambda c: [c.out("amrfinder", "amrfinder.tsv")],
    stdout_to_output=True,
    database=AMRFINDER_DB,
    threads=4,
)

mlst = Tool(
    name="mlst",
    summary="Multi-locus sequence types against PubMLST.",
    scope=Scope.SET,
    conda=("bioconda::mlst",),
    command=lambda c: ["mlst", *[str(a) for a in c.assemblies]],
    outputs=lambda c: [c.out("mlst", "mlst.tsv")],
    stdout_to_output=True,
)


# --- Relatedness ---------------------------------------------------

mashtree = Tool(
    name="mashtree",
    summary="Alignment-free tree of the whole set, built from mash distances.",
    scope=Scope.SET,
    conda=("bioconda::mashtree",),
    command=lambda c: [
        "mashtree", "--numcpus", str(c.threads), *c.args(),
        *[str(a) for a in c.assemblies],
    ],
    params=(("--genomesize", "5000000"), ("--mindepth", "5"),
            ("--kmerlength", "21"), ("--sketch-size", "10000")),
    outputs=lambda c: [c.out("mashtree", "mashtree.newick")],
    stdout_to_output=True,
    threads=8,
)

treecluster = Tool(
    name="treecluster",
    summary="Cluster assignments cut from the mashtree.",
    scope=Scope.SET,
    conda=("bioconda::treecluster",),
    needs=("mashtree",),
    command=lambda c: [
        "TreeCluster.py", "-i", str(c.out("mashtree", "mashtree.newick")),
        "-o", str(c.out("treecluster", "treecluster.tsv")), *c.args(),
    ],
    # --threshold is mandatory; these are v2's config.yaml defaults.
    params=(("--method", "max_clade"), ("--threshold", "0.05")),
    outputs=lambda c: [c.out("treecluster", "treecluster.tsv")],
)

skani = Tool(
    name="skani",
    summary="All-against-all average nucleotide identity.",
    scope=Scope.SET,
    conda=("bioconda::skani",),
    command=lambda c: [
        "skani", "triangle", "-t", str(c.threads), "--full-matrix",
        *[str(a) for a in c.assemblies],
        "-o", str(c.out("skani", "ani.tsv")),
    ],
    outputs=lambda c: [c.out("skani", "ani.tsv")],
    threads=8,
)


# --- Pangenome and phylogeny ---------------------------------------

panaroo = Tool(
    name="panaroo",
    summary="Core and accessory gene content across the set.",
    scope=Scope.SET,
    conda=("bioconda::panaroo>=1.5",),
    needs=("bakta",),
    command=lambda c: [
        "panaroo", "--clean-mode", "strict", "-a", "core", "-t", str(c.threads),
        "-o", str(c.out("panaroo")),
        "-i", *[str(c.sample_out(s, "bakta", f"{s}.gff3")) for s in c.samples],
    ],
    outputs=lambda c: [
        c.out("panaroo", "gene_presence_absence.Rtab"),
        c.out("panaroo", "core_gene_alignment.aln"),
    ],
    threads=16,
)

snp_dists = Tool(
    name="snp-dists",
    summary="Pairwise SNP distances across the core genome.",
    scope=Scope.SET,
    conda=("bioconda::snp-dists",),
    needs=("panaroo",),
    command=lambda c: [
        "snp-dists", str(c.out("panaroo", "core_gene_alignment.aln")),
    ],
    outputs=lambda c: [c.out("snp-dists", "snp-dists.tsv")],
    stdout_to_output=True,
)

fasttree = Tool(
    name="fasttree",
    summary="Approximate maximum-likelihood tree from the core genome alignment.",
    scope=Scope.SET,
    conda=("bioconda::fasttree",),
    needs=("panaroo",),
    command=lambda c: [
        "FastTree", "-nt", "-gtr",
        str(c.out("panaroo", "core_gene_alignment.aln")),
    ],
    outputs=lambda c: [c.out("fasttree", "fasttree.newick")],
    stdout_to_output=True,
    threads=4,
)


# --- Metabolic models ----------------------------------------------

carveme = Tool(
    name="carveme",
    summary="Draft genome-scale metabolic model for each genome.",
    scope=Scope.GENOME,
    conda=("bioconda::carveme",),
    needs=("bakta",),
    # Solves with open-source SCIP (pyscipopt/pulp) — no CPLEX licence needed.
    # Was commented out of v2's default target; verify it still runs, and check
    # whether DIAMOND has to be supplied explicitly (v2 shipped a separate env).
    command=lambda c: [
        "carve", str(c.out("bakta", f"{c.sample}.faa")),
        "--output", str(c.out("carveme", f"{c.sample}.xml")),
    ],
    outputs=lambda c: [c.out("carveme", f"{c.sample}.xml")],
)


CATALOGUE = Registry([
    seqkit, checkm2, gtdbtk, sylph,
    bakta, amrfinder, mlst,
    mashtree, treecluster, skani,
    panaroo, snp_dists, fasttree,
    carveme,
])

__all__ = ["CATALOGUE", "Context", "Database", "Scope", "Tool"]
