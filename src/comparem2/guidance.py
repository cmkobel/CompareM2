"""What each tool does, and how to read what it produced.

The report is the product, and a table of numbers nobody can interpret is not a
result. This module carries the explanatory half: for every tool, what question
it answers, how it works, what the specific columns on screen mean, and what
would make a conclusion drawn from them indefensible.

**Why this is not in `catalogue.py`.** A tool's spec there is ~20 lines of
executable detail — the command line, which DESIGN.md calls the largest single
risk in the rewrite. Burying that under prose written for a different audience
would hide the thing developers actually edit. The two change for different
reasons, so they live apart; `test_every_tool_has_guidance` keeps the split
honest by failing when a tool is added without its explanation.

**Sourcing rule.** Every number here was copied from the tool's own paper and
checked back against the PDF text. Where a statement is ordinary methodological
caution rather than something a paper establishes, it says so in the sentence —
a reader cannot audit a claim whose provenance is hidden, and an unattributed
caveat that sounds like a finding is worse than no caveat.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Citation:
    """One paper. `note` carries anything that qualifies what to cite."""

    text: str
    doi: str
    note: str = ""

    @property
    def url(self) -> str:
        return f"https://doi.org/{self.doi}"


@dataclass(frozen=True)
class Guidance:
    """The explanatory content for one tool.

    `reading` is the load-bearing field: pairs of (what is on screen, what it
    means). It is keyed to the columns the report actually renders rather than
    to everything the tool can emit, because guidance for output the reader
    cannot see is just more text to skip.
    """

    blurb: str  # what question this answers, for a non-bioinformatician
    method: str  # how it actually works, from the paper
    reading: tuple[tuple[str, str], ...]
    caveats: tuple[str, ...]
    citation: Citation
    # Underlying methods that deserve credit in their own right — the algorithm
    # a wrapper wraps, the database a classifier classifies against.
    also: tuple[Citation, ...] = field(default_factory=tuple)


# --- Papers ---------------------------------------------------------
# Named so a citation used by more than one tool is written once.

SEQKIT = Citation("Shen W, Le S, Li Y, Hu F (2016) SeqKit: a cross-platform and "
                  "ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE 11:e0163962",
                  "10.1371/journal.pone.0163962")
SEQKIT2 = Citation("Shen W, Sipos B, Zhao L (2024) SeqKit2: a Swiss army knife for "
                   "sequence and alignment processing. iMeta 3:e191",
                   "10.1002/imt2.191")
CHECKM2 = Citation("Chklovski A, Parks DH, Woodcroft BJ, Tyson GW (2023) CheckM2: a "
                   "rapid, scalable and accurate tool for assessing microbial genome "
                   "quality using machine learning. Nature Methods 20:1203–1212",
                   "10.1038/s41592-023-01940-w",
                   note="The PDF shipped alongside this report is the bioRxiv preprint "
                        "(doi:10.1101/2022.07.11.499243); cite the published paper.")
DIAMOND = Citation("Buchfink B, Reuter K, Drost H-G (2021) Sensitive protein alignments "
                   "at tree-of-life scale using DIAMOND. Nature Methods 18:366–368",
                   "10.1038/s41592-021-01101-x")
GTDBTK = Citation("Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH (2022) GTDB-Tk v2: "
                  "memory friendly classification with the Genome Taxonomy Database. "
                  "Bioinformatics 38:5315–5316",
                  "10.1093/bioinformatics/btac672")
GTDB = Citation("Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil P-A, Hugenholtz P "
                "(2022) GTDB: an ongoing census of bacterial and archaeal diversity "
                "through a phylogenetically consistent, rank normalized and complete "
                "genome-based taxonomy. Nucleic Acids Research 50:D785–D794",
                "10.1093/nar/gkab776",
                note="Cite alongside GTDB-Tk, and record the release — this run used r226.")
BAKTA = Citation("Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A "
                 "(2021) Bakta: rapid and standardized annotation of bacterial genomes "
                 "via alignment-free sequence identification. Microbial Genomics 7:000685",
                 "10.1099/mgen.0.000685")
PRODIGAL = Citation("Hyatt D, Chen G-L, LoCascio PF, Land ML, Larimer FW, Hauser LJ "
                    "(2010) Prodigal: prokaryotic gene recognition and translation "
                    "initiation site identification. BMC Bioinformatics 11:119",
                    "10.1186/1471-2105-11-119")
PYRODIGAL = Citation("Larralde M (2022) Pyrodigal: Python bindings and interface to "
                     "Prodigal. Journal of Open Source Software 7:4296",
                     "10.21105/joss.04296")
AMRFINDER = Citation("Feldgarden M, Brover V, Gonzalez-Escalona N, et al. (2021) "
                     "AMRFinderPlus and the Reference Gene Catalog facilitate examination "
                     "of the genomic links among antimicrobial resistance, stress "
                     "response, and virulence. Scientific Reports 11:12728",
                     "10.1038/s41598-021-91456-0")
AMRFINDER_VAL = Citation("Feldgarden M, Brover V, Haft DH, et al. (2019) Validating the "
                         "AMRFinder tool and resistance gene database by using "
                         "antimicrobial resistance genotype-phenotype correlations in a "
                         "collection of isolates. Antimicrobial Agents and Chemotherapy "
                         "63:e00483-19",
                         "10.1128/AAC.00483-19",
                         note="The source of the genotype-phenotype concordance figures.")
PUBMLST = Citation("Jolley KA, Bray JE, Maiden MCJ (2018) Open-access bacterial "
                   "population genomics: BIGSdb software, the PubMLST.org website and "
                   "their applications. Wellcome Open Research 3:124",
                   "10.12688/wellcomeopenres.14826.1",
                   note="tseemann/mlst is unpublished; this is the nomenclature it "
                        "types against, and what should be cited for an ST.")
MASHTREE = Citation("Katz LS, Griswold T, Morrison SS, et al. (2019) Mashtree: a rapid "
                    "comparison of whole genome sequence files. Journal of Open Source "
                    "Software 4:1762",
                    "10.21105/joss.01762")
MASH = Citation("Ondov BD, Treangen TJ, Melsted P, et al. (2016) Mash: fast genome and "
                "metagenome distance estimation using MinHash. Genome Biology 17:132",
                "10.1186/s13059-016-0997-x",
                note="Mashtree is a wrapper; this is where the distance measure is.")
TREECLUSTER = Citation("Balaban M, Moshiri N, Mai U, Jia X, Mirarab S (2019) TreeCluster: "
                       "clustering biological sequences using phylogenetic trees. "
                       "PLOS ONE 14:e0221068",
                       "10.1371/journal.pone.0221068")
SKANI = Citation("Shaw J, Yu YW (2023) Fast and robust metagenomic sequence comparison "
                 "through sparse chaining with skani. Nature Methods 20:1661–1665",
                 "10.1038/s41592-023-02018-3")
PANAROO = Citation("Tonkin-Hill G, MacAlasdair N, Ruis C, et al. (2020) Producing "
                   "polished prokaryotic pangenomes with the Panaroo pipeline. "
                   "Genome Biology 21:180",
                   "10.1186/s13059-020-02090-4")
CDHIT = Citation("Fu L, Niu B, Zhu Z, Wu S, Li W (2012) CD-HIT: accelerated for "
                 "clustering the next-generation sequencing data. Bioinformatics "
                 "28:3150–3152",
                 "10.1093/bioinformatics/bts565")
MAFFT = Citation("Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software "
                 "version 7: improvements in performance and usability. Molecular "
                 "Biology and Evolution 30:772–780",
                 "10.1093/molbev/mst010")
SNPDISTS = Citation("Seemann T, Klötzl F, Page AJ. snp-dists: pairwise SNP distance "
                    "matrix from a FASTA sequence alignment (software)",
                    "10.5281/zenodo.1411986",
                    note="snp-dists has no paper. There is no published description of "
                         "its algorithm, defaults or validation — cite the software DOI.")
FASTTREE = Citation("Price MN, Dehal PS, Arkin AP (2010) FastTree 2 — approximately "
                    "maximum-likelihood trees for large alignments. PLOS ONE 5:e9490",
                    "10.1371/journal.pone.0009490")
CARVEME = Citation("Machado D, Andrejev S, Tramontano M, Patil KR (2018) Fast automated "
                   "reconstruction of genome-scale metabolic models for microbial species "
                   "and communities. Nucleic Acids Research 46:7542–7553",
                   "10.1093/nar/gky537")


GUIDANCE: dict[str, Guidance] = {

    "seqkit": Guidance(
        blurb="Measures the physical shape of each assembly: how many pieces it is in, "
              "how long they are, and how much of each is G or C. It answers whether an "
              "assembly is one clean piece or a thousand fragments — nothing about "
              "whether it is complete or correct.",
        method="A single pass of `seqkit fx2tab --name --length --gc` writes one row per "
               "fasta record with its length and GC. Contig count, total length, largest "
               "contig, N50 and mean GC are all arithmetic done on those rows afterwards.",
        reading=(
            ("Contigs",
             "The number of fasta records, which for a draft assembly is a fragmentation "
             "count. A closed genome gives one record per replicon. Neither SeqKit paper "
             "offers a threshold for 'too many' — judge it against the other genomes here."),
            ("Total length",
             "Compare against the expected genome size for the taxon: far above suggests "
             "contamination or a merged bin, far below an incomplete assembly. N characters "
             "padding scaffold gaps are counted as bases, so this can overstate what was "
             "actually sequenced."),
            ("Largest / N50",
             "N50 is the length at which contigs that size or longer cover half the "
             "assembly. When Largest, N50 and Total length converge, the assembly is "
             "essentially one piece. Discarding short contigs raises N50 without improving "
             "anything, and N50 only compares between genomes of similar total size."),
            ("GC %",
             "Length-weighted, so it is the GC of the whole assembly rather than the mean "
             "of per-contig values — a small plasmid cannot drag it around. Within a "
             "species this is tight; a genome sitting well off its neighbours here is "
             "worth a second look."),
            ("The figure",
             "One bar per genome, one segment per fasta record in file order, width by "
             "length and colour by that record's own GC. Colours are relative to this run, "
             "not an absolute GC scale. A segment whose colour breaks the pattern is a "
             "prompt to check for a plasmid, phage or contaminant — not evidence of one."),
        ),
        caveats=(
            "A contiguous assembly is not a complete or uncontaminated one: completeness "
            "and contamination come from CheckM2 and identity from GTDB-Tk, not from here.",
            "The SeqKit2 paper states empty files are \"allowed and omitted without "
            "reporting errors\", so a truncated assembly yields no rows rather than a "
            "failure, and simply disappears from this table.",
            "Neither paper defines an acceptable contig count, N50 or GC deviation — they "
            "are software papers about speed and features, so any cut-off you apply comes "
            "from your own expectations for the organism.",
        ),
        citation=SEQKIT2,
        also=(SEQKIT,),
    ),

    "checkm2": Guidance(
        blurb="Answers whether a genome is good enough to use, by estimating how much of "
              "the organism's genome you recovered (completeness) and how much extra or "
              "duplicated material came with it (contamination).",
        method="Genes are called with Prodigal and annotated against a reference database "
               "with DIAMOND; those annotations plus genome length and coding statistics "
               "feed machine-learning regression models trained on simulated genomes of "
               "known quality. There is no taxonomic assignment step at any point.",
        reading=(
            ("Completeness",
             "Green is ≥90%, amber 50–90%, red <50%, matching the MIMAG bands the paper "
             "uses. Mean absolute error was 2.1±2.9% on high-quality genomes and 3.1±3.3% "
             "on medium, low-quality and highly contaminated ones — so a 92% and a 94% "
             "genome are not meaningfully different."),
            ("Contamination",
             "The predicted percentage of material that does not belong, including a second "
             "copy of something already present. Green is <5%, amber 5–10%, red ≥10%; "
             "the paper treats >10% as its own category rather than a quality tier. "
             "Mean absolute error was 1.2±1.3 on high-quality genomes."),
            ("The contamination bar",
             "Scaled ×5 by this report so that values of 1–2% stay visible. Compare the "
             "printed numbers, never the bar widths."),
            ("Green and green is not MIMAG high quality",
             "The paper is explicit that \"CheckM2 does not detect nor rely on other MIMAG "
             "factors such as full-length 16S RNA or the presence of tRNAs\". These colours "
             "cover only the completeness and contamination half of that standard."),
            ("Disagreement with a published CheckM1 number",
             "Expected, and usually CheckM2 is the one to trust. Across GTDB r202 the two "
             "agreed within 1% for 73% of completeness calls. The systematic exceptions are "
             "reduced genomes: on Patescibacteria, Dependentiae and Iainarchaeota CheckM2's "
             "error was 6.4±4.5 against CheckM1's 19.7±12.1, and on curated complete "
             "endosymbionts CheckM2 averaged 71% completeness where CheckM1 gave 39%."),
        ),
        caveats=(
            "Contamination is a predicted percentage, not a diagnosis: with no taxonomic "
            "step, CheckM2 cannot say which contigs are foreign or where they came from.",
            "Both numbers carry roughly 1–3 percentage points of residual error and the "
            "paper gives no per-genome confidence interval, so ranking genomes by small "
            "differences is not defensible.",
            "Training simulated contamination from 0% to 35%, so a value far above that "
            "range is an extrapolation — our reading of the training range, not a claim "
            "the paper makes.",
        ),
        citation=CHECKM2,
        also=(DIAMOND,),
    ),

    "gtdbtk": Guidance(
        blurb="Gives each assembly a seven-rank name, domain down to species, by placing "
              "it in the GTDB reference tree and comparing it to the nearest reference "
              "genome. The useful question it answers is how far down the ranks that name "
              "can actually be trusted.",
        method="Marker genes (120 bacterial, 122 archaeal) are extracted and aligned, and "
               "the genome is placed with pplacer into a backbone tree and then a "
               "class-level subtree; disagreement resolves to the lowest common ancestor. "
               "ANI to a species representative decides the species, and relative "
               "evolutionary divergence decides the ranks above it.",
        reading=(
            ("classification",
             "A rank left empty — a bare `s__` at the end — means GTDB-Tk would not commit "
             "there, usually because no reference was close enough on ANI. That is a "
             "result, not a failure: it says the genome is novel relative to this release. "
             "An alphanumeric name like `s__UBA1234 sp002345675` is a real GTDB species, "
             "not missing data; 77.0% of R06-RS202 species clusters carry such a name."),
            ("closest_genome_ani and closest_genome_af",
             "ANI is identity to the nearest reference; AF is how much of the two genomes "
             "could be aligned to compute it. Both matter — high ANI over a small AF is not "
             "evidence of the same species. Read 95–96% as the conventional species line "
             "and the high 80s to low 90s as a genus-level relative at best."),
            ("closest_genome_reference_radius",
             "The species-specific ANI cut-off of the reference that matched, not a global "
             "constant. Compare closest_genome_ani against this number rather than against "
             "95; GTDB treats a radius above 95% as an exceptional case."),
            ("red_value and classification_method",
             "classification_method says whether an ANI match or topological placement "
             "produced the name. When it is placement, red_value is the divergence of the "
             "placement node on a 0-to-1 scale. Treat it as an audit trail, not a score to "
             "threshold — neither paper publishes the per-rank RED intervals, so you cannot "
             "recompute the decision from this column."),
            ("msa_percent",
             "How much of the marker alignment this genome contributed. A low value means "
             "missing markers — fragmented, contaminated, or genuinely divergent — and a "
             "correspondingly weaker placement. No paper gives a cut-off, so read it "
             "comparatively and next to your CheckM2 completeness."),
        ),
        caveats=(
            "The species line is pragmatic, not biological: GTDB found no genetic "
            "discontinuity below 95% ANI, and 2.2% of species contain a genome at or above "
            "95% ANI to a different species — so 93–95% to the nearest reference is "
            "genuinely ambiguous and should not be reported as a confident species call.",
            "The name is pinned to a release and GTDB rewrites itself — 0.26% of genomes "
            "move species cluster between releases on average — so record the release or "
            "the result is not reproducible.",
            "A species-level hit is not evidence of a cultured relative: over 50% of "
            "bacterial taxa at every rank consist only of MAGs and/or SAGs, so the matched "
            "reference may itself be a bin with no phenotype behind it.",
            "The only two genuinely conflicting calls in the v1-versus-v2 comparison were "
            "both poor-quality genomes, so read a surprising name next to CheckM2 and "
            "msa_percent before trusting it.",
        ),
        citation=GTDBTK,
        also=(GTDB,),
    ),

    "bakta": Guidance(
        blurb="Works out where the genes are in each assembly and what they are called. "
              "The counts here are the structural result; the annotation itself is what "
              "AMRFinder, Panaroo and CarveMe read later in this report.",
        method="Structural calls come from specialist tools — Prodigal for protein-coding "
               "genes (via the pyrodigal binding), tRNAscan-SE, Aragorn, Infernal against "
               "Rfam, PILER-CR. Naming is where Bakta differs from Prokka: it hashes each "
               "predicted protein and looks the hash up in a pre-built database, sending "
               "only the leftovers to DIAMOND against UniRef90 then UniRef50, so no taxon "
               "has to be specified.",
        reading=(
            ("CDS",
             "For scale, on the paper's benchmark E. coli genome Bakta called 5841 CDS "
             "against 5794 from NCBI's PGAP and 5754 from Prokka. Read this next to CheckM2 "
             "completeness and assembly size: well below expectation usually means a "
             "fragmented assembly rather than gene loss, well above means a contaminated "
             "bin. Genes under 30 codons are largely invisible — Prodigal only scores open "
             "reading frames above 90 bp."),
            ("tRNA and rRNA",
             "On the benchmark genome all four annotators compared gave equal or comparable "
             "counts, so this is not where annotators differ — it mostly reflects assembly "
             "completeness. rRNA operons are near-identical repeats that routinely collapse "
             "in short-read assemblies, so a low rRNA count is more often an assembly "
             "artefact than biology (general caution, not a claim of the paper)."),
            ("Other features",
             "A sum over unlike categories — tmRNA, ncRNA genes and regulatory regions, "
             "CRISPR arrays, gaps, origins — normally dominated by non-coding RNA; the "
             "benchmark genome had 223 ncRNA genes and 66 regulatory regions. Because it "
             "mixes categories, a difference here is not interpretable on its own."),
            ("What this table does not show",
             "The functional annotation, which is the part the paper actually benchmarks: "
             "the fraction of CDS left as 'hypothetical protein', and the COG, EC and GO "
             "assignments. Those are in the Bakta output files."),
        ),
        caveats=(
            "This run uses Bakta's light database and the paper characterises only the full "
            "one (53 GB), so its functional-annotation figures — 10.6% hypothetical "
            "proteins, 94.2% of CDS identified — cannot be assumed to hold here. Since "
            "Bakta finds small proteins only by matching known sequences and \"is not able "
            "to predict these small protein coding genes de novo either\", a smaller "
            "database can only find fewer of them.",
            "CDS counts are annotator-specific: on the same 35 genomes DFAST called 127,053, "
            "Prokka 130,360 and Bakta 130,683. Compare within this run, not against a number "
            "produced elsewhere.",
            "Annotation quality collapses for organisms poorly represented in public "
            "databases — across 362 GenBank genomes of undefined genera the per-genome "
            "identification rate ran from 0 to 99.9% with a median of 10.4%, so for a MAG "
            "from an uncharacterised lineage a high CDS count says little about how much of "
            "it is functionally described.",
            "Only the assembly is passed, with no replicon metadata, so on a closed genome "
            "a gene crossing the sequence start appears as two CDS entries rather than one.",
        ),
        citation=BAKTA,
        also=(PRODIGAL, PYRODIGAL, DIAMOND),
    ),

    "amrfinder": Guidance(
        blurb="Reports which known antimicrobial-resistance genes — plus, because `--plus` "
              "is passed, biocide, metal and virulence genes — are present in each genome, "
              "against NCBI's curated Reference Gene Catalog. It tells you what resistance "
              "machinery a genome encodes, not what the isolate will do in a susceptibility "
              "test.",
        method="Predicted proteins are searched with BLASTP against curated reference "
               "proteins and with HMMER against family HMMs carrying manually curated "
               "cutoffs; a hierarchy of gene families then gives the most specific name the "
               "sequence actually supports rather than just the closest hit. This pipeline "
               "runs protein-only with `--plus` and no `--organism`.",
        reading=(
            ("Total",
             "Reported hits, not distinct phenotypes. Genes of one operon count separately — "
             "a single vancomycin-resistance cluster appears as five rows (vanR-A, vanS-A, "
             "vanH-A, vanA, vanX-A)."),
            ("Class columns",
             "The broad drug class the element contributes to; slash-joined labels mean the "
             "reference gene is annotated to all of them. Curators use a broad label "
             "deliberately where \"the literature is unclear, contradictory as to resistance "
             "phenotype, or the effect of the element is highly dependent on strain or "
             "species background\" — so breadth signals uncertainty about phenotype, not "
             "breadth of resistance."),
            ("Non-antibiotic columns",
             "Because `--plus` is on, this table mixes antibiotic resistance with stress "
             "response: the catalog held 52 biocide, 8 heat and 148 metal resistance genes. "
             "A COPPER/SILVER hit says nothing about antibiotics."),
            ("Empty cells",
             "Weaker evidence than a hit. Only genes in the catalog can be found, and "
             "carrying nothing is an ordinary result — 34.2% of the 6,242 validation "
             "isolates were pansusceptible."),
            ("Method, %Coverage and %Identity",
             "In the per-genome TSV rather than this table, and they decide how much to "
             "trust a cell. EXACTP or ALLELEP is a full-length identical match; BLASTP is "
             "above cutoff but not identical; PARTIALP means under 90% of the reference "
             "length. A PARTIALP hit at 51% coverage is a fragment just above the 50% floor, "
             "not a demonstrated functional gene."),
        ),
        caveats=(
            "A hit is a gene, not a phenotype. The authors state the mission is to identify "
            "proteins with the capacity to contribute to resistance, and that the phenotype "
            "also depends on expression and on factors outside the tool's scope such as "
            "porin mutations and efflux overexpression.",
            "The 98.4% genotype-phenotype concordance was measured across 87,679 "
            "susceptibility tests in only Salmonella, Campylobacter and E. coli, on a "
            "collection partly selected for resistant isolates — the authors say this "
            "\"might overestimate the overall PPV while underestimating the NPV\" and might "
            "not hold in other species. Per isolate rather than per test, 17% had at least "
            "one mismatch.",
            "Protein-only mode with no --organism screens no point mutations, so resistances "
            "that are mutational — Campylobacter fluoroquinolone and macrolide resistance is "
            "almost entirely gyrA T86I and 23S A2075G — are invisible here. It also cannot "
            "flag internal stops or frameshifts.",
            "Calls are made on Bakta's predicted proteins, so a gene annotation missed or "
            "truncated cannot be reported; the authors note assembly-based systems \"can be "
            "only as good as the genomic data they are assessing\".",
        ),
        citation=AMRFINDER,
        also=(AMRFINDER_VAL,),
    ),

    "mlst": Guidance(
        blurb="Puts a short, globally shared name on each genome — a sequence type, which "
              "labels the exact combination of alleles at a handful of housekeeping genes. "
              "It answers whether this genotype has been seen and named before, not how "
              "closely related two genomes are.",
        method="Each assembly is BLAST-searched against PubMLST's catalogue of known allele "
               "sequences; every locus gets the identifier of the allele it matches, and "
               "that combination is looked up in the scheme's profile table. The ST exists "
               "only if a curator has already registered that profile — these identifiers "
               "are curated nomenclature, not a computed measurement.",
        reading=(
            ("Scheme",
             "Which PubMLST scheme was auto-detected, per genome — no scheme is forced here, "
             "so different rows may carry different schemes. Two rows are only comparable "
             "if this column reads the same. PubMLST hosted schemes for over 100 species or "
             "genera, so a genome outside that coverage shows `-` and nothing useful "
             "elsewhere in the row."),
            ("Sequence type",
             "An index into that scheme's profile table, nothing more. ST 5 in one scheme "
             "has no relationship to ST 5 in another, and the numbers carry no ordering — "
             "ST 78 is not between ST 77 and ST 79. A `-` means the allele combination was "
             "not in the table: either a locus failed to call, or the combination is novel "
             "and would need submitting to PubMLST to get a number."),
            ("Alleles",
             "One entry per locus. Classical MLST indexes \"multiple, but few (six or "
             "seven), housekeeping gene fragments\", so expect roughly that many. The number "
             "in brackets is a catalogue identifier assigned at curation, not a similarity "
             "score — allele 4 is no more similar to allele 5 than to allele 100."),
            ("A dash at a locus",
             "More often a contig break than real gene loss. The paper is explicit that a "
             "core gene can be absent from a WGS dataset through \"technical issues due to "
             "incomplete\" assembly, so check contiguity before concluding a genome is "
             "untypeable. Any missing or inexact locus normally means no ST."),
        ),
        caveats=(
            "An ST is meaningless without its scheme. A scheme is an arbitrary grouping of "
            "loci chosen for a purpose, with no limit on how many exist, so STs must never "
            "be compared across schemes, species, or databases — Salmonella and E. coli are "
            "hosted on Enterobase and only mirrored into PubMLST.",
            "Six or seven housekeeping genes is a coarse ruler: the paper's own argument is "
            "that outbreak-level resolution needs cgMLST or wgMLST at ~1500–2000 loci, where "
            "point-source outbreaks are called at five or fewer allele differences. Two "
            "genomes sharing an ST are the same lineage, not the same isolate — use the ANI, "
            "SNP-distance or tree sections for relatedness.",
            "Allele and ST numbers come from a static PubMLST copy bundled with the "
            "installed tool while PubMLST itself grows continuously, so a profile with no ST "
            "today may acquire one later. Record the tool and database version alongside any "
            "published ST. (From how the tool is packaged, not from the paper.)",
        ),
        citation=PUBMLST,
    ),

    "mashtree": Guidance(
        blurb="Draws a tree of all the genomes at once without aligning anything, by "
              "compressing each assembly into a fingerprint of hashed k-mers and turning "
              "the fraction two genomes share into a distance. It answers roughly which "
              "genomes group together, in seconds.",
        method="Mash hashes every canonical 21-mer and keeps the 10,000 smallest hashes as "
               "each genome's sketch; the fraction two sketches share estimates the Jaccard "
               "index, which becomes the Mash distance D = −(1/k)·ln(2j/(1+j)). Mashtree "
               "feeds that all-against-all matrix to QuickTree's neighbour-joining.",
        reading=(
            ("Horizontal branch length",
             "Mash distance, which tracks 1 − ANI with a root-mean-square error of 0.00274 "
             "against alignment-based ANI at k=21 on 500 Escherichia genomes. As a rule of "
             "thumb from the same paper, a total path of about 0.05 corresponds to roughly "
             "95% ANI. Vertical spacing only keeps labels apart and means nothing."),
            ("No scale bar",
             "The phylogram is scaled to fit the page. Branch lengths are comparable to each "
             "other but no absolute distance can be read off the picture — the numbers are "
             "in the newick file."),
            ("No support values",
             "Not a rendering omission. Mashtree implements bootstrapping and jackknifing, "
             "but this pipeline runs plain `mashtree`, so nothing was resampled. Treat every "
             "split as untested: a group that looks tight here has not been shown to be "
             "stable."),
            ("The leftmost node",
             "A drawing convention, not an inferred ancestor — neighbour-joining output is "
             "unrooted. Do not read 'earliest-branching' off the left edge."),
        ),
        caveats=(
            "This is not a phylogeny and the authors say so: Katz et al. write \"Although "
            "Mashtree does not infer phylogeny\", and Ondov et al. \"emphasize that Mash is "
            "not explicitly designed for phylogeny reconstruction\". Use it to triage and "
            "group; use the core-genome tree for any evolutionary claim.",
            "Mash distance measures resemblance of whole k-mer sets, so it moves with gene "
            "content and genome size as well as with point mutations — Mash deliberately "
            "penalises size differences. A more fragmented or contaminated genome will shift "
            "for reasons that are not divergence, so cross-check CheckM2 and the assembly "
            "sizes before believing an odd placement.",
            "The demonstrated accuracy is for closely related genomes: the ANI correlation "
            "was shown over 90–100% ANI and degrades beyond it as the variance of the "
            "estimate grows. Deep splits in a mixed-genus set are the least trustworthy part "
            "of the drawing.",
        ),
        citation=MASHTREE,
        also=(MASH,),
    ),

    "treecluster": Guidance(
        blurb="Chops the mashtree into groups so each genome gets a cluster label instead "
              "of you having to eyeball the tree. Two genomes share a cluster only if they "
              "sit within 0.05 of each other along that tree.",
        method="It solves a min-cut partitioning problem: cut the fewest edges so every "
               "resulting group stays under a diversity limit. This pipeline uses "
               "`--method max_clade --threshold 0.05`, which limits the longest leaf-to-leaf "
               "path inside a group and additionally requires each group to be a clade — a "
               "node plus all its descendants.",
        reading=(
            ("ClusterNumber",
             "An arbitrary integer with no ordering or meaning beyond grouping. Each cluster "
             "is a clade of the mashtree whose widest internal leaf-to-leaf distance is at "
             "most 0.05. That clade requirement is what max_clade buys over plain `max`, "
             "which allows any connected piece of the tree."),
            ("ClusterNumber = −1",
             "Not an error code. TreeCluster writes −1 for every genome that ended up alone, "
             "so all −1 rows are separate singletons rather than one large cluster. "
             "Singletons are routine — 27% of clusters were singletons in the paper's 16S "
             "benchmark at threshold 0.09."),
            ("The coloured leaves in the mashtree section",
             "The same assignments, painted onto the tree they were cut from. They summarise "
             "that drawing; they are not independent evidence about it. A leaf left in the "
             "default colour means its name did not match a row here — a naming mismatch, "
             "not a biological result."),
        ),
        caveats=(
            "The 0.05 threshold is in the branch-length units of the tree it was handed — "
            "mash distances — so it means no two genomes in a cluster are more than 0.05 "
            "apart along that tree. It is not 95% ANI and not a species boundary; the paper "
            "uses a similarity threshold precisely \"to avoid the notoriously difficult "
            "problem of defining species for microbial organisms\".",
            "Clusters are an artefact of the threshold, not discovered structure: the "
            "paper's own sweep from 0.005 to 0.15 moved the answer from 181,574 to 10,112 "
            "clusters. Re-run at neighbouring thresholds and see which pairs stay together "
            "before reporting a grouping — the clustering itself takes seconds.",
            "The optimal clustering is not unique — the number of equally optimal partitions "
            "can be exponential in the number of leaves — so a genome sitting at the "
            "threshold boundary has one of several equally valid assignments.",
            "TreeCluster never looks at the genomes, only the tree, so any error in the "
            "sketching or the neighbour-joining topology is inherited whole. These clusters "
            "cannot be more reliable than the mashtree they came from.",
        ),
        citation=TREECLUSTER,
    ),

    "skani": Guidance(
        blurb="Computes average nucleotide identity between every pair of genomes — the "
              "percent identity over the parts that share ancestry. The shaded matrix is "
              "the fastest way to spot duplicates, see which inputs are the same species, "
              "and find the odd one out.",
        method="A sparse subset of k-mers is chained to locate orthologous regions, the "
               "query is cut into 20-kb chunks, and identity is estimated per chunk from the "
               "fraction of seeds that anchor into a chain, then averaged and debiased "
               "against a MUMmer-based reference. Because identity is measured only inside "
               "regions that chain, sequence missing from an incomplete assembly does not "
               "drag ANI down the way it does for pure sketching.",
        reading=(
            ("Matrix cells",
             "Estimated percent identity over the regions that could be matched. The matrix "
             "is symmetric and the value \"does not depend on the order of the inputs\". The "
             "paper refers to \"the standard 95% ANI species threshold\" as the convention "
             "for calling two genomes the same species, and describes the estimator as "
             "accurate at ANI ≥ ~82% — so this is a within-species and near-species "
             "instrument."),
            ("The shading range",
             "Stretched to the lowest value actually present, not a fixed scale. If every "
             "genome is one species the scale may run 98–100% and tiny differences will look "
             "dramatic; add one distant genome and everything else washes out. Read the "
             "numbers; colour only ranks them within this run."),
            ("Empty or near-zero cells",
             "skani declines to report a pair when the marker screen puts putative ANI below "
             "80% or the aligned fraction is under ~15%. So a blank means 'too little "
             "detectable homology to estimate', not '0% identical' — and it carries no "
             "information about how distant the pair actually is."),
            ("Differences inside the 98–100% band",
             "Against a BLAST-style baseline on 4,350 E. coli genomes skani gave Pearson R "
             "and mean absolute error of (0.981, 0.131), so gaps of one or two tenths of an "
             "ANI point are inside the method's own error. On a diverse all-to-all set "
             "spanning lower identities the same statistics were (0.976, 1.279)."),
            ("Aligned fraction, which is not shown here",
             "`skani triangle --full-matrix` writes an ANI matrix only, so this table gives "
             "identity without saying how much of each genome it covers. Since an ANI is "
             "emitted with AF as low as ~15%, a high value between a full chromosome and a "
             "small plasmid or partial MAG is possible and means something quite different. "
             "Check the assembly sizes in the QC section first."),
        ),
        caveats=(
            "ANI without aligned fraction is half the picture, and this report shows no AF "
            "column: a high identity between genomes of very different size or completeness "
            "can reflect a shared plasmid, prophage or conserved core rather than "
            "whole-genome relatedness.",
            "This runs at defaults (k = 15, c = 125), which the paper describes as tuned for "
            "similar genomes. For fragmented assemblies the authors recommend c = 70 at ANI "
            "≤95 or N50 ≤10 kb, and c = 30 near N50 ~3 kb, and report that lowering c "
            "improves both ANI and aligned-fraction accuracy — so a MAG-heavy set is being "
            "run at the least accurate of the described settings.",
            "The 95% species figure is a convention the paper cites rather than validates, "
            "and skani's debiasing is applied only above 90% putative ANI, so a pair landing "
            "at 94.8 or 95.2 should not settle a species assignment on this table alone "
            "(methodological caution beyond what the paper claims).",
        ),
        citation=SKANI,
    ),

    "panaroo": Guidance(
        blurb="Sorts every predicted gene into clusters and reports which are in all of "
              "the genomes, which in some, and which in only one. It answers what gene "
              "content these genomes share — while actively trying to undo the annotation "
              "errors that otherwise make that answer wrong.",
        method="Genes are clustered with CD-HIT at 98% identity into a graph whose nodes are "
               "orthologue clusters and whose edges join genes that neighbour each other on "
               "a contig. That context is then used to correct annotation error: merging "
               "genes translated in different frames, re-collapsing over-split families at "
               "70% identity, deleting poorly supported degree-1 nodes, and re-searching for "
               "genes the annotator missed. This pipeline runs `--clean-mode strict` with "
               "`-a core`, which also writes the core alignment that snp-dists and FastTree "
               "consume.",
        reading=(
            ("The summary line",
             "Gene clusters is the pan genome — everything seen anywhere in this set. "
             "Distinct presence patterns is how many on/off combinations actually occur; "
             "close to 1 means near-identical gene content, close to the ceiling means "
             "presence is scattered rather than tracking a few lineages."),
            ("The presence matrix",
             "One row per genome, each vertical block one pattern, width proportional to how "
             "many clusters share it. The solid block at the far left is the core; ragged "
             "middle blocks are subsets, which often follow phylogeny or a shared plasmid or "
             "prophage; blocks at the far right are genome-specific. A row visibly emptier "
             "than the others usually means a fragmented assembly, not real gene loss."),
            ("Pangenome partitions",
             "Which table you get depends on how many genomes you ran. From 20 genomes up "
             "you get the conventional Core/Soft core/Shell/Cloud bins — the ≥99% core "
             "cutoff is the paper's convention, the rest are this report's. Below 20 you "
             "get an exact count of clusters against the number of genomes sharing them, "
             "because those bin edges are fractions of N and two of the four cannot be "
             "reached on a small set: a cluster missing from one of N genomes sits at "
             "(N−1)/N, which clears 0.95 only once N ≥ 20, and a cluster in a single genome "
             "sits at 1/N, which falls under 0.15 only once N ≥ 7. Either way, one missed "
             "gene call moves a cluster out of core, so a small core with a fat accessory "
             "is more often an annotation artefact than biology."),
            ("Genes per genome",
             "Genomes of the same species should sit within a few percent of each other. "
             "Markedly below the rest points at an incomplete assembly or missed "
             "annotations; markedly above at contamination or genuinely more plasmids. Check "
             "CheckM2 before concluding anything biological."),
            ("An unexpected shared pattern",
             "The interesting row. Check whether those genomes are also neighbours in the "
             "tree; if they are not, suspect horizontal transfer or a shared mobile element."),
        ),
        caveats=(
            "`--clean-mode strict` recursively deletes poorly supported nodes, and the paper "
            "states this \"can occasionally lead to rare plasmids being removed\" — so an "
            "absent gene here is not evidence of absence, and this run cannot show that a "
            "strain lacks a rare plasmid.",
            "The paper says Panaroo \"is not recommended for metagenomic datasets\", so a "
            "pangenome built from MAGs inherits both that warning and the MAGs' "
            "incompleteness, which surfaces as an artificially small core.",
            "Core and accessory are defined only against the genomes in this run: adding or "
            "dropping one genome rewrites every partition count, and a handful of isolates "
            "does not estimate a species pangenome.",
            "Annotation error is the thing this tool exists to fix, and the scale is worth "
            "knowing: on 413 near-clonal M. tuberculosis genomes with a maximum pairwise "
            "distance of 9 SNPs — where essentially no accessory genome should exist — other "
            "tools reported 2,584 to 3,670 spurious accessory genes, 59% of the difference "
            "attributable to genes broken across contigs.",
            "Only the ≥99% core threshold comes from the paper; the Soft core, Shell and "
            "Cloud boundaries are this report's own, and are shown only when there are "
            "enough genomes for all four to be reachable.",
        ),
        citation=PANAROO,
        also=(CDHIT, MAFFT),
    ),

    "snp-dists": Guidance(
        blurb="Counts, for every pair of genomes, how many positions differ in the core "
              "gene alignment Panaroo built — the genes shared by nearly all genomes in "
              "this run. Zero means indistinguishable across that shared core.",
        method="Panaroo aligns the core clusters and concatenates them; snp-dists walks that "
               "alignment column by column and counts, per pair, the columns where the bases "
               "differ. It runs with no options, so the result is raw uncorrected counts — "
               "no evolutionary model, no distance transformation, no normalisation by "
               "alignment length.",
        reading=(
            ("The top-left header cell",
             "Not data. snp-dists writes its own version string into the first header cell "
             "and this report shows the file as it is. The remaining headers and the first "
             "cell of each row are genome names."),
            ("Off-diagonal integers",
             "A raw count of differing alignment columns over the whole concatenated core, "
             "not a percentage and not a per-Mb rate. To compare against anything outside "
             "this run, divide by the alignment length — the same 6,000 differences mean "
             "very different things over a 0.5 Mb core and a 2 Mb core."),
            ("A zero off the diagonal",
             "The two genomes are indistinguishable across the core alignment. That is what "
             "a genome included twice looks like, and also what genuinely identical isolates "
             "look like."),
            ("Its dependence on the section above",
             "The alignment length behind these counts is set by how many genes are core in "
             "this particular sample set. The Panaroo paper shows how far that can move: on "
             "the same Klebsiella pneumoniae genomes Panaroo called 3,372 core genes where "
             "Roary called 1,800. Read these counts together with the partition table."),
        ),
        caveats=(
            "Raw counts with no evolutionary model: no correction for "
            "transition/transversion bias and none for multiple substitutions at one site, "
            "so this is a ballpark similarity measure rather than an evolutionary distance. "
            "snp-dists has no publication, so this is a methodological caveat rather than a "
            "documented limitation.",
            "Recombination is not accounted for either — a single imported tract can "
            "contribute hundreds of differences at once, so a large count can mean one "
            "transfer event rather than long independent divergence.",
            "The counts depend entirely on how large the core genome is, which depends on "
            "which genomes are in the run: adding one distant or fragmented assembly shrinks "
            "the shared core for everybody. Values are not comparable between runs with "
            "different inputs.",
            "Only the core is measured, so the entire accessory genome is invisible. Two "
            "genomes can show 0 SNPs and still differ by a plasmid or a resistance cassette "
            "— check the pangenome matrix and the AMRFinder section before calling two "
            "isolates the same.",
            "Neither paper gives a threshold for 'same strain' or 'same outbreak'. Outbreak "
            "cutoffs are species-specific and come from the epidemiological literature, not "
            "from this pipeline.",
        ),
        citation=SNPDISTS,
        also=(PANAROO,),
    ),

    "fasttree": Guidance(
        blurb="Builds a phylogenetic tree from Panaroo's core-gene alignment, so relatedness "
              "is by descent rather than by overall similarity. It stays fast on large sets "
              "by searching tree space less thoroughly than a full maximum-likelihood "
              "program.",
        method="It starts from a heuristic neighbour-joining tree, improves it with "
               "minimum-evolution subtree-pruning-regrafting, then rearranges under "
               "maximum likelihood using nearest-neighbour interchanges only — never ML SPR "
               "moves, which is why the authors call it approximately-maximum-likelihood. "
               "Rate variation is handled by the CAT approximation, picking one of 20 fixed "
               "rates per site instead of integrating over a gamma distribution.",
        reading=(
            ("Branch length",
             "Expected substitutions per site under GTR+CAT, and approximate: against "
             "branch lengths re-optimised by PhyML under gamma-4 they correlated at r = 0.90 "
             "with an average difference of 13% for lengths between 0.01 and 1.0. Good "
             "enough to see which lineages are long-branched; not a calibrated divergence "
             "estimate."),
            ("Support values, in the newick rather than the figure",
             "A number before the colon at each internal node. It is SH-like local support: "
             "the current split is compared against only its two nearest-neighbour "
             "alternatives, using the Shimodaira-Hasegawa test with 1,000 resamples. It does "
             "not mean what a bootstrap value means."),
            ("The leftmost node",
             "Not a root. FastTree stores the tree with a trifurcation and the paper states "
             "the root placement is not biologically meaningful and does not affect the "
             "likelihood. If you need direction, include an outgroup."),
            ("No scale bar and no node labels in the figure",
             "For exact branch lengths or support values, open the newick file."),
        ),
        caveats=(
            "Skipping ML SPR moves costs accuracy: the paper's own benchmark puts FastTree "
            "at 84.3% of true splits against RAxML's 88.4% on 5,000-sequence protein "
            "simulations. A branching order you intend to build an argument on should be "
            "re-checked with a program that does a full search.",
            "Support values are not bootstrap values and the authors write that they \"may "
            "be biased upwards because they do not consider all of\" the alternative "
            "topologies, so a 1.00 here is weaker evidence than a 100% bootstrap. The paper "
            "warns specifically that a well-supported node sitting inside a bush of short "
            "branches is not trustworthy.",
            "Branch lengths come from the CAT approximation, and the paper states that where "
            "accurate branch lengths are essential neither CAT nor gamma-4 is sufficient. "
            "Use them for topology and relative divergence, not as a clock.",
            "The input is a concatenated core-gene alignment, so the tree describes only "
            "shared genes and assumes one history for all of them. Recombination and "
            "horizontal transfer break that assumption and the paper addresses neither "
            "(general caution, not a FastTree finding).",
            "Every benchmark in the paper is on protein families or 16S alignments; there is "
            "none on a concatenated bacterial core-gene alignment of the kind used here, and "
            "no stated minimum number of genomes below which the approximations stop being "
            "adequate.",
        ),
        citation=FASTTREE,
        also=(PANAROO,),
    ),

    "carveme": Guidance(
        blurb="Turns each genome's predicted proteins into a genome-scale metabolic model: "
              "a machine-readable network of the reactions that organism can probably run. "
              "It answers what capabilities the annotation implies, not what the organism "
              "does in the lab.",
        method="A manually curated universal bacterial model from BiGG — 4,383 reactions and "
               "2,383 metabolites — is the starting point. Your proteins are aligned with "
               "DIAMOND against 30,814 BiGG-derived sequences, alignment scores become "
               "per-reaction confidence scores through gene-protein-reaction rules, and a "
               "mixed-integer program then 'carves': keep high-scoring reactions, drop "
               "low-scoring ones, enforce connectivity so no dead ends remain.",
        reading=(
            ("Reactions",
             "Counted from every reaction element in the SBML, so transporters, exchange and "
             "spontaneous reactions and the biomass equation are included — the paper's "
             "per-organism figures count only gene-associated reactions, so this number runs "
             "larger. For scale, across 5,587 RefSeq bacterial models the average was 1,308 "
             "reactions, the smallest 238 and the largest 2,472. A few hundred is expected "
             "for a genuinely reduced genome and suspicious for a free-living organism."),
            ("Metabolites",
             "Counted from SBML species elements, which are compartment-specific: the same "
             "compound in cytosol, periplasm and extracellular space counts three times. The "
             "universal model has \"2383 metabolites (representing 1503 unique compounds)\", "
             "so do not compare this against the paper's per-organism figure of 792, which "
             "counts unique compounds."),
            ("Genes",
             "The genes from your annotation that matched a BiGG gene and survived carving. "
             "The 5,587-model average was 691, so this is a small subset of a bacterial CDS "
             "complement — the rest of the genome has no BiGG counterpart and is invisible "
             "to the model."),
            ("Differences between genomes",
             "Annotation and database differences until shown otherwise. The paper says of "
             "its own model-size statistics that \"these numbers are indicative, as they can "
             "be biased or restricted by the quality of the gene annotation and the scope of "
             "the reaction database\". Two strains annotated differently can differ by dozens "
             "of reactions with no biological difference."),
            ("What 'not gap-filled' means for use",
             "The carving problem enforces a minimum growth rate (default 0.1 h⁻¹), so the "
             "model is built to be capable of producing biomass, and in the paper's benchmark "
             "four of five organisms grew on minimal medium from genome data alone with no "
             "gap-filling — the fifth needed one added reaction. Use these for network-level "
             "comparison and as FBA starting points; do not report a growth rate, auxotrophy "
             "or cross-feeding interaction without first checking the model grows on a "
             "defined medium."),
        ),
        caveats=(
            "The authors call these drafts themselves: models \"should still be considered "
            "as drafts subject to further refinement, as they might require "
            "organism-specific curation to reproduce certain phenotypes\".",
            "Coverage is bounded by BiGG — primary metabolism is \"essentially complete\" "
            "while \"peripheral pathways associated with secondary metabolism contain "
            "multiple gaps\" — so a missing secondary-metabolite pathway is not evidence the "
            "organism lacks it.",
            "Transport is the known weak point: the paper attributes poorer "
            "substrate-utilisation performance to \"the lack of annotated transporters\", "
            "which makes uptake, secretion and any cross-feeding claim the least trustworthy "
            "part of these models.",
            "No template flag is passed, so the generic biomass composition is used even for "
            "archaea and Gram-positives, and the archaeal template that \"contains ether "
            "lipids in the cell membrane composition and lacks peptidoglycans\" is not "
            "selected.",
        ),
        citation=CARVEME,
        also=(DIAMOND,),
    ),
}


def citations(names: list[str] | tuple[str, ...]) -> list[Citation]:
    """Every distinct paper behind `names`, primary first, in tool order.

    Deduplicated by DOI, because the underlying methods are shared — DIAMOND is
    cited by three tools and Panaroo by three.
    """
    seen: dict[str, Citation] = {}
    for name in names:
        entry = GUIDANCE.get(name)
        if entry is None:
            continue
        for citation in (entry.citation, *entry.also):
            seen.setdefault(citation.doi, citation)
    return list(seen.values())
