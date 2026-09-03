"""The v3 tool set.

Thirteen tools. Selected 2026-09-01 as fourteen; sylph was removed 2026-09-02
when reading its paper showed the command could never have worked. See DESIGN.md
for what was dropped and why.

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

import sys

from . import carve_scip, steps
from .tools import Context, Database, Registry, Scope, Tool

# --- Databases -----------------------------------------------------
# Sizes are measured, not estimated. `None` means nobody has measured it yet.

# Two tool specs are named here because their *database* fetch runs the tool's
# own binary, so the spec has to be one string used in both places. Bakta is the
# case that makes this mandatory rather than tidy: db 6.x is required by bakta
# 1.12.x and refused by nothing else, so a pin that drifted between the tool and
# its download would fetch a database the tool then rejects at runtime.
_BAKTA_SPEC = "bioconda::bakta>=1.10"
_AMRFINDER_SPEC = "bioconda::ncbi-amrfinderplus"

# What a plain URL fetch needs. Present in the pixi environment and on any
# Linux; declared because under `--use-conda` a download rule gets exactly the
# environment it asks for and nothing else.
_FETCH_TOOLS = ("conda-forge::curl", "conda-forge::tar")

_CHECKM2_URL = "https://zenodo.org/records/14897628/files/checkm2_database.tar.gz?download=1"

CHECKM2_DB = Database(
    name="checkm2",
    url=_CHECKM2_URL,
    size=1_735_095_710,  # measured 2026-09-01
    # The URL pins a Zenodo record, so the layout inside the tarball is pinned
    # with it: `CheckM2_database/uniref100.KO.1.dmnd`. Verified against an
    # existing extracted copy 2026-09-02. Linking it to a stable name means
    # catalogue.py's path does not move if a later record renames the file —
    # but if that happens, this link target has to change too.
    fetch=lambda db: [
        ["mkdir", "-p", str(db / "checkm2")],
        ["curl", "-fSL", "--retry", "3", "-o", str(db / "checkm2" / "db.tar.gz"), _CHECKM2_URL],
        ["tar", "-xzf", str(db / "checkm2" / "db.tar.gz"), "-C", str(db / "checkm2")],
        ["ln", "-sf", "CheckM2_database/uniref100.KO.1.dmnd",
         str(db / "checkm2" / "checkm2.dmnd")],
        ["rm", "-f", str(db / "checkm2" / "db.tar.gz")],
    ],
    conda=_FETCH_TOOLS,
    ready="checkm2/checkm2.dmnd",
)

# r232, because that is what the tool demands — not a preference. GTDB-Tk
# 2.7.2's `config/common.py` reads `COMPATIBLE_REF_DATA_VERSIONS = ['r232']`,
# so it refuses anything else outright. This was r226 for a day, which the
# installed tool would have rejected after a 141.4 GB download.
#
# The version coupling runs both ways and is the same shape as bakta's: pinning
# the tool without moving the database, or the database without moving the pin,
# is a runtime failure rather than a solve failure.
#
# The canonical host, and also the fastest measured from Denmark, which was not
# the expected result. Single-stream, 2026-09-02, same 60,806,405,195-byte
# object on all three:
#
#   data.gtdb.ecogenomic.org (primary)  9.7 MB/s
#   data.ace.uq.edu.au                  4.7 MB/s
#   data.gtdb.aau.ecogenomic.org        3.0 MB/s   <- the Danish mirror
#
# The cap is per connection, not per server: three extra streams alongside a
# running download each still got 3.4 MB/s off Aalborg, 13.2 MB/s in total. So
# a parallel fetcher would help — but one stream off the primary already
# reaches 9.7 of the ~14 MB/s this path can carry, which is not worth an aria2
# dependency for.
_GTDB_URL = (
    "https://data.gtdb.ecogenomic.org/releases/release232/232.0/"
    "auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
)

GTDB_DB = Database(
    name="gtdb",
    url=_GTDB_URL,
    # 60.8 GB, measured from content-length 2026-09-02 — less than half of
    # r226's 141.4 GB. GTDB-Tk 2.7 replaced FastANI and Mash with skani, so the
    # package no longer carries the reference genomes that dominated it; the
    # first directory in the tarball is now `skani/database/sketches.db`.
    size=60_806_405_195,
    # Published by the mirror at `<release>/MD5SUM.txt`. Not checked by the
    # fetch — hashing 60.8 GB takes minutes — but it is the way to tell a
    # truncated download from a corrupt one when the extraction fails.
    md5="25a59e0352b1fd150c589f56559767d4",
    # --strip-components=1 drops the tarball's own `release232/` wrapper so
    # GTDBTK_DATA_PATH can be the directory we chose rather than a name that
    # changes every release. Verified 2026-09-02 by listing the first entries
    # of the stream — `curl -r 0-3000000 <url> | tar -tzf -` — rather than by
    # downloading 60.8 GB to find out.
    #
    # The stamp remains the `ready` file: the wrapper is known now, but no
    # interior filename is, and a stamp is honest about that.
    # `-C -` because this is a six-hour transfer at the 3.0 MB/s measured from
    # thylakoid, and `--retry` alone does not resume: measured 2026-09-02, a
    # retried transfer truncated the file and restarted from byte 0. The mirror
    # sends `accept-ranges: bytes`, and a re-run of this rule now continues a
    # partial tarball instead of re-fetching 60.8 GB.
    #
    # The cost is one narrow case: if the tarball is *complete* but the
    # extraction failed, `-C -` asks for a range past the end and the mirror
    # answers 416, so curl exits 22 and the rule fails until the tarball is
    # deleted by hand. That is loud and costs one `rm`; the alternative was
    # silent and cost thirteen hours.
    fetch=lambda db: [
        ["mkdir", "-p", str(db / "gtdb")],
        ["curl", "-fSL", "--retry", "3", "-C", "-",
         "-o", str(db / "gtdb" / "db.tar.gz"), _GTDB_URL],
        ["tar", "-xzf", str(db / "gtdb" / "db.tar.gz"), "-C", str(db / "gtdb"),
         "--strip-components=1"],
        ["rm", "-f", str(db / "gtdb" / "db.tar.gz")],
        ["touch", str(db / "gtdb" / ".fetched")],
    ],
    conda=_FETCH_TOOLS,
)

# v2 used `--type full` (30 GB compressed / 84 GB on disk). v3 uses light.
# 1.3 GB is Bakta's documented figure, not a content-length measurement,
# because there is no static URL — `bakta_db download` picks the record.
#
# The database version is tied to the Bakta version: 1.12.x requires db 6.x and
# refuses 5.1 outright. So the download must run from the same environment as
# the tool, and pinning Bakta without re-fetching the database is a runtime
# failure rather than a solve failure.
BAKTA_DB = Database(
    name="bakta-light",
    url="bakta_db download --type light",
    size=None,
    # `bakta_db download` insists on creating its own `db-light` subdirectory,
    # so it is fetched into a scratch directory and moved to where the tool
    # expects it. Note the database name and the directory differ, which is why
    # `ready` is relative to the database root rather than derived from `name`.
    fetch=lambda db: [
        ["bakta_db", "download", "--output", str(db / ".bakta_dl"), "--type", "light"],
        ["rm", "-rf", str(db / "bakta")],
        ["mv", str(db / ".bakta_dl" / "db-light"), str(db / "bakta")],
        ["rm", "-rf", str(db / ".bakta_dl")],
    ],
    # `bakta_db` is Bakta's own script, so this fetch needs Bakta itself.
    conda=(_BAKTA_SPEC,),
    ready="bakta/version.json",  # verified to exist in a real db v6.0 light
)

AMRFINDER_DB = Database(
    name="amrfinder",
    url="amrfinder -u",
    size=None,
    # This one cannot go under --databases. `amrfinder -u -d <dir>` exits with
    # "AMRFinder update option (-u/--update) only operates on the default
    # database directory. The -d/--database option is not permitted" —
    # verified 2026-09-02. So the data lands in $CONDA_PREFIX and all we can
    # record here is that the update ran.
    fetch=lambda db: [
        ["amrfinder", "-u"],
        ["mkdir", "-p", str(db / "amrfinder")],
        ["touch", str(db / "amrfinder" / ".updated")],
    ],
    # The fetch *is* the tool, so under `--use-conda` the data lands in the
    # environment Snakemake built for *this rule* — and the analysis rules have
    # to end up in that same environment or they will not find it. They do,
    # because Snakemake addresses a deployed environment by
    # md5(realpath(envs_dir) + env file content) (read from conda.py in 9.26.1),
    # so byte-identical env files under one --conda-prefix are one directory on
    # disk. Sharing the spec string is what keeps them byte-identical.
    conda=(_AMRFINDER_SPEC,),
    ready="amrfinder/.updated",
    out_of_tree=True,
)


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
    # Pinned because the database is: only 2.7+ accepts r232, and the URL above
    # is r232. An older build would install cleanly and refuse the data.
    conda=("bioconda::gtdbtk>=2.7",),
    # Each genome lives in its own directory once inputs are canonicalised, so
    # --genome_dir cannot be used and GTDB-Tk takes a --batchfile instead.
    #
    # Both of the steps around this command were described in a comment here
    # for a day without existing, which is worse than not having them: the rule
    # named a batchfile nothing created, so the command would have failed on its
    # first line, and it declared a summary file the tool never writes under
    # that name, so the job would have failed on a missing output even if the
    # command had worked. Neither had ever been executed, so neither showed up.
    # No `--skip_ani_screen`: that flag does not exist in 2.7.2. It belonged to
    # the versions whose ANI screen used Mash and needed either `--mash_db` or
    # permission to skip it; 2.7 screens with skani from the reference package
    # and the flag is gone. Checked against the installed tool's own `--help`,
    # which is the only place this was discoverable — the drafted command would
    # have died on `unrecognized arguments`.
    command=lambda c: [
        "gtdbtk", "classify_wf", "--cpus", str(c.threads),
        "--batchfile", str(c.workdir / ".comparem2" / "gtdbtk_batchfile.tsv"),
        "--out_dir", str(c.out("gtdbtk")),
    ],
    # Two columns, tab separated, no header: FASTA path, then the genome id.
    # The id is the sample name, which is already wildcard-safe, and which is
    # what makes the summary joinable to every other tool's output.
    #
    # Kept in `.comparem2/` beside the Snakefile rather than in the tool's own
    # output directory: it is a generated file of the same kind, and putting a
    # rule's input inside the directory the tool writes to invites the tool to
    # remove it.
    files=lambda c: {
        c.workdir / ".comparem2" / "gtdbtk_batchfile.tsv": "".join(
            f"{c.sample_out(s, f'{s}.fna')}\t{s}\n" for s in c.samples),
    },
    # GTDB-Tk writes `<prefix>.bac120.summary.tsv` and `<prefix>.ar53.summary.tsv`
    # separately, and an all-bacterial set produces no ar53 file at all. Both
    # candidate directories are passed because the documentation does not say
    # which one classify_wf uses, and it has differed between versions.
    # The step is named by absolute script path, not `-m comparem2.steps`: see
    # steps.py. `-m` needs the package importable, and pixi provides it through
    # a relative PYTHONPATH that does not survive a rule's working directory.
    post=lambda c: [
        [sys.executable, steps.__file__, "merge-tsv",
         "--out", str(c.out("gtdbtk", "gtdbtk.summary.tsv")),
         str(c.out("gtdbtk", "*.bac120.summary.tsv")),
         str(c.out("gtdbtk", "*.ar53.summary.tsv")),
         str(c.out("gtdbtk", "classify", "*.bac120.summary.tsv")),
         str(c.out("gtdbtk", "classify", "*.ar53.summary.tsv"))],
    ],
    outputs=lambda c: [c.out("gtdbtk", "gtdbtk.summary.tsv")],
    # GTDB-Tk takes its database location *only* from the environment; there is
    # no flag for it. Without this the `--databases` value was silently ignored
    # and the tool would have used whatever GTDBTK_DATA_PATH happened to hold,
    # or failed. Found 2026-09-02 while wiring up the downloads.
    env=lambda c: [("GTDBTK_DATA_PATH", str(c.databases / "gtdb"))],
    database=GTDB_DB,
    threads=16,
)

# sylph was here, and is gone — see DESIGN.md 2026-09-02. It is a read-based
# metagenome profiler: `sylph profile` reads samples from FASTQ or .sylsp and
# treats positional FASTA as *reference genomes*, so passing assemblies to it
# could never have worked. The deeper problem is that v3's input is assemblies,
# which is not the question sylph answers; skani already covers fast
# assembly-to-assembly identity.


# --- Annotation and screening --------------------------------------

bakta = Tool(
    name="bakta",
    summary="Structural and functional genome annotation.",
    scope=Scope.GENOME,
    conda=(_BAKTA_SPEC,),
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
    conda=(_AMRFINDER_SPEC,),
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
        *c.args(),
        *[str(a) for a in c.assemblies],
        "-o", str(c.out("skani", "ani.tsv")),
    ],
    # skani's default c = 125 is tuned for complete, similar genomes. Its paper
    # documents three presets: c = 200 for >95% ANI with N50 >10 kb, c = 70 for
    # ANI <=95 or N50 <=10 kb, and c = 30 for N50 ~3 kb. v3 is a wide view that
    # has to survive fragmented MAGs, and it cannot know which it was given, so
    # it takes the middle setting: lowering c costs runtime and index size but
    # improves both ANI and aligned-fraction accuracy, and skani is fast enough
    # (>50x FastANI) to absorb it. Override with `--set skani-c=125` for a set
    # of complete isolate genomes.
    params=(("-c", "70"),),
    # `triangle --full-matrix` writes two files: the ANI matrix at -o, and the
    # aligned-fraction matrix at that path plus `.af`. Declaring both is not
    # bookkeeping — an ANI is emitted once alignment covers as little as ~15%
    # of a genome, so identity without coverage is half the result, and a high
    # ANI between a chromosome and a partial MAG means something quite
    # different. Verified from a real run's log 2026-09-02.
    outputs=lambda c: [c.out("skani", "ani.tsv"), c.out("skani", "ani.tsv.af")],
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
    # One thread, and measured rather than assumed. `bioconda::fasttree` ships
    # both `FastTree` and `FastTreeMP`, and `Tool.env` could hand the latter its
    # OMP_NUM_THREADS — so the switch is reachable. It just does not pay: on a
    # real 7-taxon x 2,066,459 bp core alignment, plain took 242.3 s and 209.4 s
    # on two reps, and FastTreeMP took 251.7 / 278.5 / 254.6 / 204.3 / 208.2 s
    # at 1 / 2 / 4 / 8 / 16 threads — every one inside the spread of the plain
    # binary against itself. On the 4-taxon test set MP is 2.3% slower.
    # Nor does that change with more genomes, which is the obvious objection:
    # 8 threads bought 1.13x at a simulated 25 taxa, 1.15x at 100 and 1.09x at
    # 400, for 40-65% more CPU and 8 cores Snakemake would reserve. Numbers in
    # DECISIONS.md, 2026-09-03. Re-measure before changing this.
    threads=1,
)


# --- Metabolic models ----------------------------------------------

carveme = Tool(
    name="carveme",
    summary="Draft genome-scale metabolic model for each genome.",
    scope=Scope.GENOME,
    # DIAMOND, pyscipopt and scip come with it — checked against the bioconda
    # 1.6.6 recipe, which matters because `carve` shells out to DIAMOND and an
    # environment holding only carveme would fail at its first step.
    conda=("bioconda::carveme",),
    needs=("bakta",),
    # Solves with open-source SCIP, so no CPLEX licence — but *not* with the
    # presolver conda-forge's SCIP ships. Measured 2026-09-02: 601 s in the MILP
    # and a model missing 253 of its annotated reactions, against 43 s and 45
    # with `presolving/milp/maxrounds=0`. Neither CarveMe nor ReFramed exposes a
    # solver parameter, so it is set by a wrapper; carve_scip.py carries the
    # numbers and the evidence that the presolver is wrong here rather than
    # slow.
    #
    # Run by a bare `python`, deliberately. The wrapper imports carveme, so it
    # needs the interpreter of the environment *the tool* is in — which under
    # `--use-conda` is not the one running CompareM2. That is the mirror image
    # of the `sys.executable` in gtdbtk's post step, which runs our own code.
    command=lambda c: [
        "python", carve_scip.__file__,
        "--faa", str(c.out("bakta", f"{c.sample}.faa")),
        "--output", str(c.out("carveme", f"{c.sample}.xml")),
        *c.args(),
    ],
    # Two, because `carve` names its DIAMOND output after its input and the
    # wrapper gives it an input inside this directory. Undeclared, that file
    # landed in Bakta's directory instead and silently replaced Bakta's own
    # `<sample>.tsv` — see carve_scip.py.
    outputs=lambda c: [
        c.out("carveme", f"{c.sample}.xml"),
        c.out("carveme", f"{c.sample}.tsv"),
    ],
    # argv[0] is an interpreter now, and the preflight would check `python`,
    # find it everywhere, and never report carveme missing.
    executable="carve",
)


CATALOGUE = Registry([
    seqkit, checkm2, gtdbtk,
    bakta, amrfinder, mlst,
    mashtree, treecluster, skani,
    panaroo, snp_dists, fasttree,
    carveme,
])

__all__ = ["CATALOGUE", "Context", "Database", "Scope", "Tool"]
