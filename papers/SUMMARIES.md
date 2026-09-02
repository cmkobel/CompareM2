# Paper summaries for the v3 tool set

Extracted 2026-09-02 by reading each tool's paper(s) in full alongside the exact
command line in `src/comparem2/catalogue.py` and the columns rendered in
`src/comparem2/report.py`.

**Every number below was copied verbatim from its paper and verified by
substring match against `pdftotext` output — 178 of 178 passed.** The condensed
version of this file is `src/comparem2/guidance.py`, which is what the report
renders. This file is the long form, and it keeps two things the report drops:
the full quantitative record, and the *Not established* notes recording what
could not be determined from the papers. Those notes are the honest limit of
what the guidance can claim.

Ordering follows the `CATALOGUE` registry.

---

## seqkit
*Shen W, Sipos B, Zhao L 2024, iMeta 3:e191, doi:10.1002/imt2.191 (SeqKit2); earlier version: Shen W, Le S, Li Y, Hu F 2016, PLoS ONE 11(10):e0163962, doi:10.1371/journal.pone.0163962*

SeqKit measures the basic physical shape of each assembly: how many pieces (contigs) it is
in, how long each one is, and how much of each is G or C. It answers "is this genome one
clean piece or a thousand fragments, and does its GC content look like a single organism?" —
nothing more.

**How it works.** CompareM2 runs a single SeqKit subcommand, `seqkit fx2tab --name --length --gc`, which the
2016 paper lists as "Converting FASTA/Q to tabular format with extra information": it walks
the FASTA once and writes one row per record with the record name, its length in bases, and
its GC content. Every figure in the report's table — contig count, total length, largest
contig, N50, mean GC — is arithmetic that CompareM2 does on those rows afterwards; SeqKit
itself only supplies the per-contig lengths and GC.

### Reading the output

**Contigs.** The number of FASTA records in the file, nothing else. A closed genome gives one record per
replicon (chromosome plus one per plasmid); a draft assembly gives one per contig, so this
column is really a fragmentation count. Neither SeqKit paper offers a threshold for 'too
many' — judge it against the other genomes in the same run and against how the assembly was
produced.

**Total length.** The sum of all record lengths, in bases. Compare it with the expected genome size for the
taxon (which has to come from elsewhere in the report, e.g. the taxonomy or completeness
sections) — a total far above expectation suggests contamination or a duplicated bin, far
below suggests an incomplete assembly. Note that N characters used to pad scaffold gaps are
counted as bases here.

**Largest / N50.** Largest is the longest single record. N50 is computed by sorting records longest-first and
taking the length at which the running sum first reaches half of Total length, i.e. half the
assembly sits in contigs at least that long. When Largest, N50 and Total length are all
close, the assembly is essentially one piece. Read N50 together with Contigs and Total
length: discarding short contigs raises N50 without improving anything, and N50 is only
comparable between genomes of similar total size.

**GC %.** The length-weighted mean of the per-record GC values, so it is the GC of the whole assembly
rather than the average of the individual contig values (a 5 kb plasmid cannot drag it
around). Within a species this number is tight; a genome whose GC sits well off its
neighbours in the same run is worth a second look.

**Contig figure (bars coloured by GC).** One bar per genome, one segment per FASTA record in the order the records appear in the
file, segment width proportional to length, colour set by that record's own GC. The legend's
minimum and maximum are taken from all contigs across the run, so the colours are relative
to this report, not an absolute GC scale. The figure is drawn so that a fragmented assembly,
an outsized genome, or a contig whose GC does not match its neighbours stands out; treat the
last of those as a prompt to check for a plasmid, phage or contaminant, not as evidence of
one — neither SeqKit paper makes any such claim.

### Caveats

- SeqKit only describes the file it was given: a low contig count and a high N50 mean the
  assembly is contiguous, not that it is complete, correct, or derived from a single
  organism — completeness and contamination come from CheckM2, and taxonomic identity from
  GTDB-Tk, elsewhere in this report.

- The SeqKit2 paper states that empty files "are allowed and omitted without reporting
  errors", so a truncated or empty assembly yields no rows instead of a failure, and the
  report simply omits that genome from the table rather than flagging it.

- The 2016 paper states that sequence type (DNA/RNA/Protein) is detected from the leading
  subsequences of the first record only, so the GC column is meaningful only if the input
  really is a nucleotide assembly.

- Not a claim either paper makes, but relevant to reading the table: the length column
  counts every character in the record, so N characters padding scaffold gaps inflate Total
  length, Largest and N50 relative to the bases actually sequenced.

### Quantitative record (4 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `38` | SeqKit2 offers 38 subcommands; CompareM2 calls exactly one of them (fx2tab), so most of what you will read about SeqKit online does not apply to this report section. | “SeqKit2 comprises 38 subcommands” |
| `24 chromosomes, one mitochondrial sequence and 169 scaffolds` | FASTA parsing was benchmarked on a whole human genome, i.e. far larger and more fragmented than a bacterial assembly — running this step on microbial genomes is cheap. | “dataset B is the human genome with 24 chromosomes, one mitochondrial sequence and 169 scaffolds” |
| `67,748 DNA sequences with average length of 41 Kb` | The 2016 benchmarks also covered a FASTA file with tens of thousands of records, so a highly fragmented MAG is within the tested regime. | “Dataset A consists of 67,748 DNA sequences with average length of 41 Kb” |
| `780 Mb` | Peak memory plateaued as input files grew, because it is set by the length of the longest single sequence record rather than by file size (2016 tests on repeated human chromosome 1). | “780 Mb (Fig 3A and 3B)” |

### Not established from the papers (4)

- Neither paper defines N50 or offers any threshold for an acceptable contig count, N50,
  total length or GC deviation — they are software papers about speed and features, not
  assembly-quality criteria, so any cut-off you apply has to come from your own expectations
  for the organism.

- Neither paper documents how `--gc` handles ambiguous or degenerate bases (whether N is in
  the denominator, whether S and W contribute), so the exact GC definition behind the column
  is unverified here.

- Neither paper states whether `fx2tab --gc` emits a percentage (0-100) or a fraction;
  CompareM2's report renders the value directly as a percentage.

- Neither paper reports a runtime for `fx2tab` specifically — the published benchmarks cover
  sequence summary, read/write and reverse-complement tasks.

---

## checkm2
*Chklovski A, Parks DH, Woodcroft BJ, Tyson GW 2022, bioRxiv preprint (not peer reviewed at this version), doi:10.1101/2022.07.11.499243. Companion aligner: Buchfink B, Reuter K, Drost H-G 2021, Nature Methods 18:366-368, doi:10.1038/s41592-021-01101-x*

CheckM2 answers "is this genome good enough to use?" by estimating, for each assembly, how
much of the organism's genome you actually recovered (completeness) and how much extra or
duplicated material got mixed in (contamination). Both are given as percentages, one row per
genome.

**How it works.** CheckM2 calls genes with Prodigal, annotates the predicted proteins with KEGG orthology ids
using DIAMOND blastp against a bundled reference database, then feeds those annotations plus
genome length, number of coding sequences and amino acid counts into machine-learning
regression models trained on hundreds of thousands of synthetic genomes of known quality.
Completeness comes from either a gradient-boosted-tree "general" model or a neural-network
"specific" model, chosen per genome by a cosine-similarity proxy for taxonomic novelty;
contamination always comes from the gradient boost model. There is no taxonomic
classification step.

### Reading the output

**Completeness (%).** The predicted fraction of the organism's genome present in this assembly. Green in the
report is >=90%, amber 50-90%, red <50%, matching the MIMAG bands the paper uses (high
quality >90% complete, medium 50-90%, low <50%). Treat it as an estimate with real error:
benchmarking on simulated RefSeq r202 genomes gave a mean absolute error of 2.1+/-2.9% for
high-quality genomes and 3.1+/-3.3% for medium/low-quality and highly contaminated ones. A
92% and a 94% genome are not meaningfully different.

**Contamination (%).** The predicted percentage of extra material — sequence that does not belong to the target
genome, including a second copy of material already there. It is a single number, not a list
of suspect contigs, and because CheckM2 does no taxonomic assignment it cannot tell you
where the contaminant came from. Green is <5%, amber 5-10%, red >=10%; the paper treats >10%
contaminated as its own 'high contamination' category rather than a quality tier. Mean
absolute error was 1.2+/-1.3 on high-quality genomes and 1.7+/-1.7% on medium and low
quality genomes.

**Contamination bar length.** The coloured bar next to the contamination number is scaled x5 by the report so that values
of 1-2% remain visible. Read the printed number, not the bar width, when comparing genomes.

**Green + green does not equal MIMAG high quality.** MIMAG's high-quality definition also requires rRNA genes and a tRNA count. The paper is
explicit that 'CheckM2 does not detect nor rely on other MIMAG factors such as full-length
16S RNA or the presence of tRNAs' — so the colours here reflect only the
completeness/contamination half of the standard.

**Disagreement with a published CheckM1 number.** Across GTDB r202, 73% of completeness calls agreed within 1% and 91% within 5%;
contamination agreed within 1% for 82% of genomes and within 5% for 99%. The systematic
exceptions are small or reduced genomes — for lineages such as Patescibacteria, Dependentiae
and Iainarchaeota CheckM2's completeness error was 6.4+/-4.5 against CheckM1's 19.7+/-12.1,
and on curated complete endosymbionts CheckM2 averaged 71% completeness against CheckM1's
39%. If your genome is a small-genome symbiont or a CPR/DPANN lineage, expect CheckM2 to
read higher, and the paper argues it is the more accurate of the two.

### Caveats

- A green completeness and a green contamination value do not make a genome MIMAG high
  quality: the paper states plainly that "CheckM2 does not detect nor rely on other MIMAG
  factors such as full-length 16S RNA or the presence of tRNAs", so the rRNA and tRNA
  requirements of the standard are simply not assessed here.

- Contamination is a predicted percentage, not a diagnosis — CheckM2 runs no taxonomic
  assignment step, so it cannot say which contigs are foreign or what organism they came
  from; the paper notes it may somewhat underestimate contamination from divergent (higher-
  taxon) sources in medium-quality genomes, while being far less likely than CheckM1 to
  overestimate.

- Both numbers are model predictions with residual error of roughly 1-3 percentage points
  (completeness MAE 2.1±2.9% high quality, 3.1±3.3% medium/low; contamination MAE 1.2±1.3
  high quality), so ranking genomes by small differences is not defensible; the paper
  reports cohort averages and gives no per-genome confidence interval.

- The models were trained on complete bacterial and archaeal RefSeq genomes with simulated
  self-contamination spanning 0% to 35%, so a reported contamination far above that range
  sits outside the trained interval — that extrapolation concern is our reading, not a
  statement the paper makes.

### Quantitative record (16 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `>90% complete, <5% contaminated` | MIMAG high-quality band used to bin all benchmarking results | “(high quality: >90% complete, <5% contaminated; medium quality:” |
| `medium 50-90% complete <10% contaminated; low <50% complete <10% contaminated` | MIMAG medium- and low-quality bands | “50%-90% complete, <10% contaminated; low quality: <50% complete, <10% contaminated)” |
| `>10% contaminated` | Genomes above 10% contamination are treated as a separate category, not a quality tier | “a separate group for high contamination (>10% contaminated)” |
| `MAE 2.1±2.9%` | Completeness error on high-quality simulated RefSeq r202 genomes | “CheckM2 and CheckM1 on high quality genomes (CheckM2 MAE: 2.1±2.9%, CheckM1 MAE:” |
| `MAE 3.1±3.3% (CheckM1 4.7±5.4%)` | Completeness error on medium, low-quality and highly contaminated genomes | “genomes (Figure 3a; CheckM2 MAE: 3.1±3.3%, CheckM1 MAE: 4.7±5.4%)” |
| `MAE 1.2±1.3` | Contamination error on high-quality genomes | “mean average error (MAE: 1.2±1.3) was comparable” |
| `1.7±1.7% (CheckM1 3.0±4.0%)` | Contamination error on medium and low quality genomes | “medium and low quality genomes (CheckM2: 1.7±1.7%, CheckM1: 3.0±4.0%)” |
| `73% within 1% of each other` | Agreement with CheckM1 completeness across GTDB r202 genomes | “with 73% of all completeness predictions being within 1% of each other,” |
| `82% within 1%, 99% within 5%` | Agreement with CheckM1 contamination across GTDB r202 genomes | “for contamination, with 82% of all genome predictions being within 1% of each other, and 99%” |
| `MAE 6.4±4.5` | Completeness error on reduced-genome lineages (Patescibacteria, Dependentiae, Iainarchaeota) | “completeness predictions are far more accurate (MAE: 6.4±4.5) than” |
| `MAE 19.7±12.1` | CheckM1 completeness error on the same reduced-genome lineages | “CheckM1 (MAE: 19.7±12.1; Figure 4b)” |
| `71% (CheckM1: 39%)` | Average completeness predicted for curated complete bacterial endosymbiont genomes | “completeness, with CheckM2 predicting an average completeness of 71%, compared to” |
| `0% to 35%` | Contamination range spanned by the self-contamination training simulation | “to generate a range of 0% to 35% contamination” |
| `5% to 100% at 5% intervals` | Completeness range spanned by the training simulation | “completeness between 5% to 100% at 5% intervals” |
| `50% completeness` | GTDB's minimum completeness cutoff, which the paper argues is inaccurate when derived from CheckM1 | “minimum cutoff (50% completeness)” |
| `1.56±0.83 genomes per minute per thread` | Throughput benchmarked against CheckM1 | “1.56±0.83 genomes per minute per thread on an AMD EPYC 7702 64-Core Processor” |

### Not established from the papers (6)

- The preprint gives no per-genome uncertainty estimate or quality flag, so nothing in the
  rendered table tells you how confident CheckM2 is about that particular row.

- No upper limit is stated above which contamination predictions should be distrusted; the
  cross-contamination benchmark capped the contaminating fraction at 10%, and training self-
  contamination reached 35%.

- The paper describes the cosine-similarity rule that selects the 'general' (gradient boost)
  versus 'specific' (neural network) completeness model per genome, but I did not verify
  whether CompareM2's quality_report.tsv column for that choice is surfaced anywhere in the
  report — the report renders only Genome, Completeness and Contamination.

- I read only the bioRxiv preprint (CC-BY-ND, explicitly not peer reviewed at that version).
  Whether any of these numbers changed in the peer-reviewed Nature Methods article is
  unverified here.

- The throughput figure was measured with 45 threads on an AMD EPYC 7702; CompareM2 runs
  `checkm2 predict --threads 8` by default, so expected wall time on your machine cannot be
  read off it directly.

- Behaviour on inputs outside the training scope — eukaryotic bins, viral contigs, plasmid-
  only assemblies — is not addressed by the paper.

---

## gtdbtk
*Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH 2022, Bioinformatics 38(23):5315-5316, doi:10.1093/bioinformatics/btac672 — with the underlying taxonomy described in Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil P-A, Hugenholtz P 2022, Nucleic Acids Research 50(D1):D785-D794, doi:10.1093/nar/gkab776*

GTDB-Tk gives each of your assemblies a full seven-rank name — domain down to species — by
placing it in the Genome Taxonomy Database reference tree and comparing it to the closest
reference genome. It answers "what is this genome, and how far down the ranks can that name
actually be trusted?"

**How it works.** Marker genes are extracted and aligned (GTDB uses 120 bacterial and 122 archaeal markers),
and the genome is placed with pplacer first into a backbone tree holding one representative
per family, then re-placed into the appropriate class-level subtree; if the two disagree,
"the genome is classified by taking the lowest common ancestor between the backbone and
class-level subtree". The final rank calls "use the same RED and ANI criterion as GTDB-Tk
v1": average nucleotide identity to a species representative decides the species, and
relative evolutionary divergence — how far down the reference tree the placement node sits —
decides the ranks above it.

### Reading the output

**classification.** The full string d__;p__;c__;o__;f__;g__;s__. A rank left empty (e.g. a bare 's__' at the
end) means GTDB-Tk would not commit at that rank — usually because no reference was close
enough on ANI, so the call fell back to tree position. That is a result, not a failure: it
says your genome is novel relative to this GTDB release. Equally, an alphanumeric name like
s__UBA1234 sp002345675 is a genuine GTDB species, not a placeholder for missing data — 77.0%
of GTDB R06-RS202 species clusters carried such a name. Note the report prints only the
first 50 genomes of this table.

**closest_genome_ani and closest_genome_af.** ANI is the average nucleotide identity to the nearest reference genome; AF is the alignment
fraction, i.e. how much of the two genomes could be aligned to compute it. Both must clear a
bar for a species call — a high ANI over a small AF is not evidence of the same species.
GTDB built its own clusters with an AF criterion of 0.65, which the paper says it expected
to relax to 0.5 to accommodate incomplete MAGs. Read a value of 95-96% as the conventional
species line ('An ANI threshold of 95-96% is widely considered the gold standard'), and
anything in the high 80s to low 90s as a genus-level relative at best.

**closest_genome_reference_radius.** The species-specific ANI cut-off of the reference your genome matched, not a global
constant. GTDB's analysis treated 95% as the normal value and noted that species with a
radius above 95% 'are exceptional cases'. Compare closest_genome_ani against this number,
not against 95, when deciding whether the species call was earned.

**red_value, classification_method and note.** classification_method tells you which criterion produced the name — an ANI match to a
reference, or topological placement. When it is placement, red_value carries the relative
evolutionary divergence of the placement node, a 0-to-1 scale where the root is 0 and the
tips are 1, used to decide which rank the node corresponds to. Treat red_value as an audit
trail rather than a score to threshold: neither paper publishes the per-rank RED intervals
the classifier applies, so you cannot recompute the rank decision from this column alone.

**msa_percent and warnings.** msa_percent is how much of the marker-gene alignment your genome contributed (of 120
bacterial or 122 archaeal markers in the GTDB R06-RS202 tree). A low value means missing
markers — fragmented assembly, heavy contamination, or a genome that is genuinely divergent
— and the placement it feeds is correspondingly weaker. Neither paper states a numeric cut-
off, so use it comparatively across your own genomes and read it next to your CheckM2
completeness. Any text in 'warnings' should be read before the classification it sits
beside.

### Caveats

- The species line is pragmatic, not biological: GTDB's own analysis found no genetic
  discontinuity below 95% ANI ("nearly equal numbers of pairs with ANI values between 78%
  and 95%"), and 2.2% of species contain a genome at or above 95% ANI to a different species
  — so a genome landing at 93-95% ANI to its nearest reference is genuinely ambiguous and
  should not be reported as a confident species call.

- The name is pinned to a GTDB release (CompareM2 downloads r226), and GTDB rewrites itself:
  on average 96.7% of species representatives are unchanged between releases and 0.26% of
  genomes move to a different species cluster, so always record the release alongside the
  classification or the result is not reproducible.

- Classification is only as good as the assembly: the sole two genuinely conflicting calls
  in the GTDB-Tk v1-vs-v2 comparison were "both relatively poor-quality genomes", so read
  this table next to CheckM2 completeness/contamination and the msa_percent column before
  trusting a surprising name.

- A species-level hit is not evidence of a cultured relative — over 50% of bacterial taxa at
  every rank consist exclusively of MAGs and/or SAGs, so the reference your genome matched
  may itself be a metagenome-assembled bin with no phenotype behind it.

### Quantitative record (22 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `320 GB` | GTDB-Tk v1 needed this much RAM to place genomes into the R07-RS207 bacterial reference tree with pplacer | “320 GB of RAM” |
| `<55 GB` | The divide-and-conquer placement in v2 cut peak memory requirements | “320 GB to <55” |
| `12 genomes (0.07%)` | Of 16 710 GEM bacterial genomes representing novel taxa, this many were classified differently by v2 than v1 | “Only 12 genomes (0.07%) did not” |
| `13 genomes (0.06%)` | In a second evaluation on 23 548 dereplicated genomes, this many differed between v1 and v2 | “Only 13 genomes (0.06%) had different GTDB-Tk v1 and GTDB-Tk v2 classifications” |
| `six over-classified, four under-classified` | The few incongruent calls were mostly off by one rank in either direction | “over(six genomes) or under-classified (four genomes) by a single” |
| `2 genomes` | The only two genomes with genuinely conflicting v1/v2 assignments were poor-quality inputs | “these were both relatively poor-quality genomes” |
| `22-35%` | Speed-up of v2 over v1 when processing 1000 genomes on 1-64 CPUs | “v2 also ran 22–35% faster when processing 1000 genomes with” |
| `95-96%` | Conventional ANI range used to delineate prokaryotic species | “An ANI threshold of 95–96% is widely considered” |
| `78% and 95%` | GTDB found no genetic discontinuity below the species line: intra-genus representative pairs are spread evenly across a wide ANI range instead of clustering | “we find nearly equal numbers of pairs with ANI values between 78% and 95%” |
| `6.0%` | Fraction of closest between-species genome pairs that nevertheless exceed the species ANI line | “only 6.0% of the closest interspecific genome” |
| `80.8%` | Fraction of closest between-species pairs that are well separated on ANI | “pairs have an ANI ≥ 95% and 80.8% of pairs are well” |
| `2.2%` | Proportion of GTDB species that are 'fuzzy' — they contain a genome at or above 95% ANI to a genome of another species | “In terms of species, 2.2% can be described as fuzzy as” |
| `89.5%` | Proportion of GTDB species that are cleanly discrete (no genome above 94% ANI to another species) | “89.5% can be described as discrete as they” |
| `77.0%` | Proportion of GTDB R06-RS202 species clusters whose name is an alphanumeric placeholder rather than a Latin binomial | “NCBI lacking a species classification and 77.0% of GTDB” |
| `120 bacterial and 122 archaeal` | Marker genes used to infer the domain-specific GTDB reference trees (R06-RS202) | “The 120 bacterial and 122 archaeal” |
| `0.65` | Alignment-fraction criterion used to assign genomes to GTDB species clusters, which GTDB expected to relax for incomplete MAGs | “change the AF criteria used for assigning genomes to GTDB species clusters from 0.65 to” |
| `96.7%` | Average proportion of GTDB species representatives that stay the same from one release to the next | “RS04-RS89 with 96.7% of representatives being unchanged” |
| `0.26%` | Average proportion of genomes that move to a different species cluster between GTDB releases | “Only 0.26% of genomes on average are assigned to a different species cluster between releases” |
| `254 090 bacterial and 4316 archaeal genomes` | Size of the GTDB R06-RS202 reference set | “254 090 bacterial and 4316 archaeal genomes” |
| `Over 50%` | Share of bacterial taxa at any rank that consist only of MAGs and/or SAGs, i.e. have no cultured representative | “Over 50% of bacterial taxa, regardless of rank, consist exclusively of MAGs” |
| `95%` | The typical ANI circumscription radius of a GTDB species cluster | “ANI circumscription radius of 95%” |
| `FastANI v1.3` | Tool GTDB used to compute the ANI and AF values that define its species clusters | “ANI and AF values were calculated with FastANI v1.3” |

### Not established from the papers (6)

- Neither paper mentions an ANI pre-screen or the --skip_ani_screen flag that CompareM2
  passes — both predate GTDB-Tk 2.2, where the screen was introduced. I could not establish
  from these papers what skipping it costs. The intended reading (a shortcut that lets
  genomes matching a known species representative bypass marker alignment and pplacer
  placement, so skipping it trades runtime and memory for nothing else) is not something
  either paper supports, and should be verified against the GTDB-Tk documentation before
  being printed as fact.

- Neither paper publishes the numeric RED intervals GTDB-Tk uses to convert a placement node
  into a rank, so red_value in the output table cannot be interpreted against a stated
  threshold.

- Neither paper states a minimum msa_percent, a minimum AF, or a minimum ANI that GTDB-Tk
  itself applies at classification time. The 0.65 AF figure quoted above is for constructing
  GTDB species clusters, not for GTDB-Tk classification, and should not be presented as a
  GTDB-Tk cut-off.

- Neither paper covers GTDB r226, the release CompareM2 downloads. The figures here are from
  R06-RS202 and R07-RS207. In particular the archaeal marker set has changed since the 122
  markers reported by Parks 2022 — the pipeline's own code comments refer to ar53 output
  files — so the marker count for the run at hand may differ.

- The exact column set of gtdbtk.summary.tsv is version-dependent and I had no example
  output on disk to check against; the columns named in the guidance above are the standard
  GTDB-Tk v2 summary fields and should be confirmed against a real run.

- The papers do not describe what GTDB-Tk writes when it declines a rank (empty field vs.
  absent rank prefix), so the operational reading of a truncated classification string is
  inferred, not cited.

---

## sylph
*Shaw J & Yu YW 2024, Nature Biotechnology 43:1348–1359, doi:10.1038/s41587-024-02412-y*

Sylph asks, for every genome in a reference collection, "how much of this genome's DNA is
present in my sample, and at what nucleotide identity?" — and it answers in seconds rather
than hours, which is why CompareM2 runs it as a provisional taxonomic call while GTDB-Tk is
still working. It reports a match only when the identity is above the 95% species boundary,
so a row here is a species-level statement about the nearest thing already in the database.

**How it works.** Sylph subsamples k-mers (k = 31, roughly one in 200) from each reference genome and from
your sample, then measures what fraction of a reference's k-mers is contained in the sample;
raising that containment to the power 1/k gives the "naive" containment ANI. Its actual
contribution is a correction: it models the k-mer counts as a zero-inflated Poisson, infers
the effective coverage λ, and divides the containment by (1 − e^−λ) before taking the k-th
root, so that k-mers missing merely because the sample did not cover the genome are not
mistaken for sequence divergence; `profile` additionally assigns each shared k-mer to
whichever genome has the highest ANI (winner-take-all), recomputes ANI, and prints only
genomes above 95%.

### Reading the output

**Genome_file.** The reference genome sylph matched, printed as its filename or accession in the GTDB sketch
— not a species name. CompareM2's command does not run sylph's optional taxonomy-mapping
step, so you have to look the accession up. Note also that `profile` assigns each shared
k-mer to the single highest-ANI genome, so among several near-identical genomes of the same
species the particular one listed is close to a tie-break: trust the species, not the
strain.

**Adjusted_ANI vs Naive_ANI.** Naive_ANI is the raw containment raised to 1/31; Adjusted_ANI is the same quantity after
dividing containment by (1 − e^−λ), i.e. what the identity would have been had the sample
covered the genome fully. Rows only appear at all above 95% adjusted ANI. The gap between
the two columns is the size of the coverage correction: on 90%-identity Nanopore reads the
naive estimate fell below 95% while sylph 'still gave good estimates (>99% median ANI)'. For
a complete assembly the two should sit close together; a large gap means the correction, not
the data, is doing the work.

**Eff_cov, Eff_lambda, Median_cov, Mean_cov_geq1.** Effective coverage λ is the depth of the reference's k-mers actually observed after read
errors and read length are accounted for, so it is lower than true sequencing depth. It is
what makes low-abundance detection possible: for a genome whose true containment ANI was
99.26%, sylph recovered >95% ANI from as little as 0.008x effective coverage. The λ
correction is only applied when the median k-mer multiplicity is ≤3; above that, the
reported coverage is a robust mean or the median itself.

**ANI_5-95_percentile and Lambda_5-95_percentile.** A 90% interval built by resampling the k-mer multiplicities over 100 iterations and taking
the 5th and 95th percentiles. A wide ANI interval means λ was poorly determined and the
point estimate is soft. The authors caution that these intervals are optimistic: 'confidence
intervals are slightly tight but the coverage probabilities are not less than 70% under
simulations.'

**Taxonomic_abundance and Sequence_abundance.** Taxonomic abundance is each detected genome's λ divided by the sum of all detected λ (a per-
organism share); sequence abundance weights each λ by genome length (a per-base share), so
long genomes take a larger share of it than of the taxonomic column. For a clean single-
isolate assembly you expect one row near 100 in both. Extra rows carrying real abundance
mean k-mer content from a second species is present — check that against CheckM2
contamination before calling it a mixed culture; the paper does not test sylph on assemblies
for this purpose.

### Caveats

- The command CompareM2 currently runs (`sylph profile <db>.syldb <assembly>.fna ... -o
  profile.tsv`) passes the assemblies as positional FASTA arguments, which `sylph profile`
  treats as additional reference genomes rather than as samples; with no FASTQ or .sylsp
  sample present sylph exits with 'No read files found' and this section will be empty —
  verified by reading sylph's source (src/contain.rs, src/cmdline.rs), not by running it,
  and not a point the paper makes.

- Every result in the paper is for shotgun reads, and the zero-inflated Poisson correction
  is a model of read sampling; an assembly carries each k-mer essentially once, so the
  coverage columns do not mean what they mean in the paper and Adjusted_ANI should be read
  alongside Naive_ANI rather than instead of it (this is our reading of the paper's stated
  filters, not something the paper tests).

- The paper states plainly that sylph 'loses sensitivity when organisms only have database
  representatives at higher taxonomic ranks', so a genome that is a novel species can
  produce no sylph row at all while GTDB-Tk still returns a genus or family — silence from
  sylph is not evidence of a bad assembly.

- Containment ANI is not the quantity GTDB-Tk compares against its species representatives;
  the paper notes it is biased slightly upward at k = 31 with the bias subsiding only above
  95% ANI, so small numeric disagreement between sylph's ANI and GTDB-Tk's is expected
  rather than a sign that one of them is wrong.

### Quantitative record (14 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `>95%` | Sylph reports only genomes above the species-level containment ANI cutoff | “outputs genomes for which ANI is >95% by default.” |
| `k = 31` | k-mer length used for sketching | “Sylph first subsamples the k-mers (k = 31) for each genome in a” |
| `c = 200` | k-mer subsampling rate (roughly one k-mer in c is kept) | “c = 200 by default, effectively speeding up computation and lowering” |
| `>50` | Minimum number of sampled k-mers before a genome is reported at all | “of FracMinHash k-mers for a genome is >50 by default, although this” |
| `0.008x` | Lowest effective coverage at which the correction still recovered a >95% ANI call, in the downsampled K. pneumoniae test | “However, beginning at even 0.008×” |
| `99.26%` | Reference value the correction was aiming at in that test (ANI from the full, undownsampled read set) | “the true containment ANI (99.26% as estimated with all reads)” |
| `median k-mer multiplicity ≤3` | Condition under which the coverage adjustment is applied at all | “for ANI if the median multiplicity for k-mers X1, …, XN is ≤3” |
| `100 iterations` | Bootstrap iterations behind the ANI percentile column | “over 100 iterations and performs coverage-adjusted ANI calculations” |
| `5th and 95th percentile` | Percentiles reported as the ANI interval | “Sylph then takes the 5th and 95th percentile ANI” |
| `not less than 70%` | How well calibrated those intervals are in simulation | “confidence intervals are slightly tight but the coverage probabilities are not less than 70% under simulations.” |
| `>99% median ANI` | Median adjusted ANI on error-prone (90% identity) Nanopore reads, where naive ANI dropped below the 95% cutoff | “sylph still gave good estimates (>99% median ANI).” |
| `92% mean precision, 82% F1` | Species-level precision and F1 on a synthetic community where most organisms had no species-level database representative | “sylph had 92% mean precision and 82% F1” |
| `>95%` | ANI above which containment ANI stops being noticeably biased upward relative to standard genome-to-genome ANI at k = 31 | “this bias subsides for high ANI (>95%)” |
| `<4 GB of memory for >25,000 genomes` | Memory footprint and database size in the profiling benchmark, the reason sylph can finish while GTDB-Tk is still running | “time while taking <4 GB of memory for >25,000 genomes.” |

### Not established from the papers (6)

- Which GTDB release the sylph sketch is built from is still undecided in the pipeline
  (SYLPH_DB is marked url="TBD" in catalogue.py), so a name-level disagreement with GTDB-Tk
  caused purely by different reference releases cannot be ruled out.

- Whether Adjusted_ANI saturates at 100 when the sample is an assembly rather than reads:
  the paper's filters (λ adjustment requires ≥3 k-mers in each of the two most-populated
  multiplicity classes) suggest the adjustment is usually skipped for assemblies, but sylph
  is not installed in this environment and no run was made to check.

- The paper gives no threshold for when a second sylph row in one sample should be called
  contamination, and no guidance at all on using assembled genomes as sylph samples.

- The output column names quoted here were read from sylph's source (print_header in
  src/contain.rs); the paper does not list them, and they were not confirmed against a real
  profile.tsv.

- The paper contains no benchmark of sylph against GTDB-Tk, so nothing here supports a claim
  about how often the two agree on isolate genomes.

- CompareM2's report currently renders sylph with the generic fallback table, which shows at
  most the first 50 lines of profile.tsv (header plus 49 rows) with a 'First 50 rows.' note;
  with many input genomes the later samples will simply not appear.

---

## bakta
*Schwengers et al. 2021, Microbial Genomics 7:000685, doi:10.1099/mgen.0.000685*

Bakta reads each assembly and works out where the genes are and what they are called. The
table shows how many protein-coding genes (CDS), tRNAs, rRNAs and other features Bakta found
in each genome; the full annotation it writes (GFF3 plus a protein FASTA) is what AMRFinder,
Panaroo and CarveMe read later in this report.

**How it works.** Structural annotation is done by separate specialist tools: protein-coding genes are called
de novo by Prodigal (Hyatt 2010, run through the pyrodigal binding, Larralde 2022), tRNAs by
tRNAscan-SE, tmRNA by Aragorn, rRNA and other non-coding RNAs by Infernal against Rfam
covariance models, and CRISPR arrays by PILER-CR. Naming the proteins is where Bakta differs
from Prokka: instead of aligning everything, it takes an MD5 hash of each predicted protein
sequence and looks it up in a pre-built database of known proteins ("alignment-free sequence
identification"), sending only the leftovers to DIAMOND against UniRef90 clusters and then
UniRef50 — no taxon needs to be specified, because the database is not taxon-specific.

### Reading the output

**CDS.** Protein-coding genes called in that genome. For scale, on the paper's benchmark genome (E.
coli O26:H11 strain 11368) Bakta called 5841 CDS, against 5794 from NCBI's PGAP, 5754 from
Prokka and 5740 from DFAST. Read this column next to the CheckM2 completeness and the
assembly size: a count well below what the taxon should carry usually means a fragmented or
partial assembly rather than gene loss, and a count well above usually means a contaminated
bin. Genes shorter than 30 codons are largely invisible here — Prodigal only scores start-
stop pairs above 90 bp.

**tRNA.** tRNAscan-SE calls. The paper reports that on its benchmark genome all four annotators (PGAP,
Prokka, DFAST, Bakta) predicted equal or comparable numbers of tRNAs, tmRNAs, rRNAs and
CRISPR arrays, so this column is not where annotators differ — it mostly reflects how
complete and contiguous your assembly is. Neither paper gives an expected tRNA count for a
complete genome, so treat an unusually low value as a prompt to check the assembly, not as a
measured deficiency.

**rRNA.** Ribosomal RNA genes found by Infernal using Rfam covariance models. Bakta reported 22 rRNA
features on the benchmark E. coli genome, the same number as PGAP, Prokka and DFAST. General
caution, not a claim the paper makes: rRNA operons are near-identical repeats and routinely
collapse in short-read assemblies, so a low rRNA count is more often an assembly artefact
than biology.

**Other features.** Everything else Bakta writes to the GFF3 except the one contig record per sequence — tmRNA,
ncRNA genes, ncRNA regulatory regions, CRISPR arrays, gaps, origins of replication and
origins of transfer. It is normally dominated by non-coding RNA: on the benchmark genome
Bakta called 223 ncRNA genes and 66 regulatory regions, plus 2 CRISPR arrays and 4 origins
of replication. Because it is a sum over unlike categories, a difference between genomes in
this column is not interpretable on its own — open the GFF3 and count by feature type.

**What the table does not show.** The functional annotation, which is the part of Bakta the paper actually benchmarks: the
fraction of CDS left as 'hypothetical protein', and the COG/EC/GO terms and RefSeq/UniRef
cross-references. Those are in the Bakta output files, not in this table. In this pipeline
the protein FASTA feeds AMRFinder and CarveMe and the GFF3 feeds Panaroo, so the CDS count
here sets the input size for those sections.

### Caveats

- CompareM2 v3 downloads Bakta's light database, and this paper characterises only the full
  one (53 GB, database v3.0) — none of its functional-annotation figures (10.6 %
  hypothetical proteins, 94.2 % of CDS identified) can be assumed to hold for the annotation
  in front of you, and because the paper states plainly that Bakta "is not able to predict
  these small protein coding genes de novo either" but finds sORFs only by matching known
  sequences, a smaller database can only find fewer small proteins.

- CDS counts are annotator-specific, so compare them within this run and not against a
  number produced elsewhere: on the same 35 genomes DFAST called 127 053, Prokka 130 360 and
  Bakta 130 683 CDS, and the paper attributes most of the Bakta-Prokka gap to 235 detected
  small proteins plus differences in the internal feature overlap filters.

- Annotation quality drops sharply for organisms poorly represented in public databases — on
  362 GenBank genomes of undefined genera the per-genome identification rate ranged from 0
  to 99.9 % with a median of 10.4 %, so for a MAG from an uncharacterised lineage a high CDS
  count says little about how much of it is functionally described.

- CompareM2 passes Bakta only the assembly — no replicon metadata on completeness or
  topology. The paper describes that metadata as the mechanism by which Bakta merges partial
  CDS pairs running off the 5' and 3' edges of a complete replicon, so on a closed genome a
  gene crossing the sequence start will appear as two CDS entries rather than one.

### Quantitative record (19 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `5841 (Bakta) vs 5794 (PGAP)` | CDS called by Bakta on the benchmark genome E. coli O26:H11 strain 11368, versus PGAP | “Bakta (n=5841) and PGAP (n=5794) predicted more genes” |
| `223` | ncRNA genes called by Bakta on the benchmark genome — the bulk of the 'Other features' column | “highest number of ncRNA genes (n=223)” |
| `4` | Origins of replication called by Bakta on the benchmark genome; no other tool called any | “predicting origins of replication (n=4)” |
| `82 (Bakta), 44 (PGAP)` | Small proteins (sORFs) detected on the benchmark genome by Bakta and by PGAP; neither Prodigal nor MetaGeneAnnotator predicts these de novo | “(n=82) and PGAP (n=44) that are not predicted de novo by” |
| `90 bp` | Minimum open reading frame length Prodigal will score, which sets the practical floor on the CDS column | “every start-stop pair above 90 bp in the entire genome.” |
| `10.6 %` | CDS left as 'hypothetical protein' by Bakta across 35 taxonomically diverse RefSeq genomes (full database) | “annotated as hypothetical protein as low as 10.6 % (n=13 902)” |
| `31.6 % and 41.2 %` | Same 35 genomes, hypothetical-protein ratio for DFAST and Prokka | “which achieved total ratios of 31.6 and 41.2 %, respectively” |
| `94.2 %` | Predicted CDS identified exactly against known protein sequences by alignment-free lookup, across the 35 RefSeq genomes | “able to identify publicly known UPSs; 94.2 % (n=123 105) of” |
| `24.2 %` | Hypothetical-protein ratio on 198 high-quality bacterial MAGs (CheckM completeness >=95, contamination <=1) | “as 24.2 % (n=138 282) outperforming DFAST (n=232 516)” |
| `38.6 %` | CDS identified by alignment-free lookup on those 198 MAGs — far lower than for RefSeq isolates | “Bakta was able to precisely identify 38.6 % (n=220 753) of all” |
| `0 to 99.9 %, median 10.4 %` | Genome-wise identification rate on 362 GenBank genomes of undefined genera — the worst case for a poorly represented organism | “between 0 and 99.9 %, respectively, with a median of 10.4 %.” |
| `25.4 %` | Hypothetical-protein ratio on those 362 genomes of undefined genera | “protein as low as 25.4 % (n=286 406), outperforming DFAST” |
| `0.9 %` | How much the alignment-free lookup contributes to functional annotation quality, versus running DIAMOND alignments alone — it buys speed and cross-references, not better products | “as low as 0.9 % (n=1164)” |
| `at least 80 % and 90 %` | DIAMOND filters for assigning a protein to a UniRef90 cluster: mutual coverage and sequence identity | “of at least 80 and 90 %, respectively” |
| `7:09 min:s` | Wall-clock runtime for one E. coli genome on eight threads (CompareM2 also runs Bakta with 8 threads) | “were considerably shorter than those of Bakta (7:09 min:s).” |
| `4.4 GB` | Peak memory for one genome | “memory than Bakta (4.4 GB). Also, database sizes of Prokka” |
| `53 GB` | Size of the Bakta database the paper benchmarks (v3.0, the full build) — CompareM2 v3 uses the light build instead | “that of Bakta (53 GB)” |
| `1 to 16 cores` | Range over which Bakta gains from extra CPU cores | “tool exhibited a solid scalability between 1 and 16 CPU cores.” |
| `identical predictions` | pyrodigal, the Prodigal binding Bakta calls, versus the Prodigal executable | “resulting in identical predictions” |

### Not established from the papers (5)

- The paper (2021, Bakta 1.1, database v3.0) predates the light/full database split
  entirely. It does not say what the light build contains, how large it is, or how much
  functional annotation, dbxref assignment or sORF detection is given up relative to full —
  I could not source any of that from the papers provided.

- Neither paper gives an expected or acceptable range for tRNA or rRNA counts per genome, so
  there is no published threshold behind a judgement that a value in those columns is 'too
  low'.

- Whether Bakta writes detected sORFs into the GFF3 as CDS features (and therefore whether
  the 82 small proteins in the benchmark sit inside the CDS column or the Other-features
  column of this table) is not stated in the paper, and I did not check a real output file.

- The benchmarks are bacterial throughout — the MAG set was explicitly screened for "a
  taxonomical assignment within the bacterial GTDB lineage" — so the paper provides no
  evidence about annotation quality on archaea, which CompareM2 also accepts.

- Neither paper reports how much the CDS count itself (as opposed to the functional
  descriptions) changes with database size, so I cannot say whether the light database
  shifts the numbers in this table at all.

---

## amrfinder
*Feldgarden M, Brover V, Gonzalez-Escalona N, et al. 2021, Scientific Reports 11:12728, doi:10.1038/s41598-021-91456-0; validation paper: Feldgarden M, Brover V, Haft DH, et al. 2019, Antimicrobial Agents and Chemotherapy 63:e00483-19, doi:10.1128/AAC.00483-19*

AMRFinderPlus takes the proteins Bakta predicted in each genome and reports which known
antimicrobial-resistance genes — plus, because CompareM2 passes --plus, known
biocide/metal/acid stress-resistance and virulence genes — are present, against NCBI's
manually curated Reference Gene Catalog. It tells you what resistance machinery a genome
encodes; it does not tell you what the isolate will do in a susceptibility test.

**How it works.** Every predicted protein is searched with BLASTP against the curated reference proteins and
with HMMER against family HMMs that have manually curated cutoffs; matches below the per-
gene curated identity cutoff (default drops < 90% identity) or below the coverage cutoff
(default drops < 50% of the reference length) are discarded, and a hierarchy of gene
families is used to give the most specific name the sequence actually supports instead of
simply naming the closest hit. CompareM2 runs protein-only mode with --plus and no
--organism, so stress and virulence families are included but point mutations, internal-
stop/frameshift detection and taxon-specific suppression of ubiquitous genes are all off —
the paper calls the combined nucleotide+protein mode "the most accurate method".

### Reading the output

**Total.** The number of reported hits for that genome, not the number of distinct resistance
phenotypes. Genes of one operon are counted separately: in a CompareM2 test run on
Enterococcus faecium the five GLYCOPEPTIDE hits in genome E8202 were vanR-A, vanS-A, vanH-A,
vanA and vanX-A — one vancomycin-resistance cluster, five rows.

**Class columns (MACROLIDE, GLYCOPEPTIDE, TRIMETHOPRIM, ...).** Class is the broad drug class the element contributes resistance to; slash-joined labels
such as LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN mean the reference gene is annotated to all of
them. The curators use the broad class label deliberately where 'the literature is unclear,
contradictory as to resistance phenotype, or the effect of the element is highly dependent
on strain or species background' — so a broad label signals uncertainty about phenotype, not
breadth of resistance.

**Non-antibiotic class columns (COPPER/SILVER, METAL, BIOCIDE, ...).** Because --plus is on, the table mixes antibiotic resistance with stress-response genes. The
2021 catalog held 2 acid, 52 biocide, 8 heat and 148 metal resistance genes. A COPPER/SILVER
hit (copB appears in all three genomes of the CompareM2 test set) is a metal-resistance gene
and says nothing about antibiotics.

**Empty cells (·).** A dot means no hit of that class was reported, which is weaker evidence than a hit. Only
genes in the catalog can be found, and the authors could not recover both copies of a
duplicated blaCMY-2 from draft assemblies that the closed plasmids carried. Carrying nothing
is also an ordinary result: 34.2% of the 6,242 validation isolates were pansusceptible.

**Method / % Coverage / % Identity (in amrfinder.tsv, not in the report table).** The per-genome TSV carries the evidence the table omits, and it is what decides how much to
trust a cell. EXACTP/ALLELEP is a full-length identical match; BLASTP is above cutoff but
not identical (by default > 90% identity over more than 90% of the sequence); PARTIALP means
less than 90% of the reference length inside a contig. In the CompareM2 test set, fexB in
genome 116_2 was PARTIALP at 51.81% coverage — a fragment just above the 50% floor, not a
demonstrated functional gene.

### Caveats

- A hit is a gene, not a phenotype: the 2019 authors state the tool's mission is to identify
  proteins that have the capacity to contribute to resistance, and that the actual phenotype
  depends on expression level and on factors outside its coverage such as porin mutations
  and efflux overexpression, so this table cannot replace susceptibility testing.

- The 98.4% concordance figure was measured only in Salmonella, Campylobacter and E. coli
  against CLSI/NARMS breakpoints, on a collection partly selected for resistant isolates —
  the authors say this choice 'might overestimate the overall PPV while underestimating the
  NPV' and that the agreement 'might not occur in other species', so it should not be
  carried over to environmental isolates or MAGs.

- CompareM2 runs protein-only without --organism, so no point mutations are screened; in the
  same validation set Campylobacter fluoroquinolone and macrolide resistance was almost
  entirely mutational (gyrA T86I, 23S A2075G), and those resistances are simply invisible
  here — protein-only mode also cannot flag internal stop codons or frameshifts, which is
  exactly why AMRFinder missed 16 loci that ResFinder found.

- Calls are made on Bakta's predicted proteins, so a gene that annotation missed or
  truncated cannot be reported, and the authors note assembly-based systems 'can be only as
  good as the genomic data they are assessing' — in their draft assemblies a duplicated
  blaCMY-2 copy present in the closed plasmid was lost.

### Quantitative record (13 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `98.4% of 87,679 susceptibility tests` | Overall agreement between AMRFinder's gene-based resistance predictions and measured susceptibility phenotypes, over 6,242 NARMS isolates | “Of 87,679 susceptibility tests performed, 98.4% were consistent with predictions” |
| `PPV 0.955, NPV 0.992` | Predictive values for the whole validation set: a predicted-resistant call was right 95.5% of the time, a predicted-susceptible call 99.2% of the time | “Overall consistency was 98.4% of susceptibility tests performed, with a PPV of 0.955 and an NPV of 0.992.” |
| `98.0% consistent, PPV 0.94, NPV 0.992` | Salmonella enterica (5,425 isolates) genotype-phenotype consistency and predictive values | “Overall consistency was 98.0%, with a PPV of 0.94 and an NPV of 0.992.” |
| `96.7% consistent, PPV 0.904, NPV 0.982` | Campylobacter coli was the worst-performing species in the validation set | “Overall consistency was 96.7%, with a PPV of 0.904 and an NPV of 0.982.” |
| `6 of 11` | Worked example of a resistance gene without the phenotype: azithromycin predictions in S. enterica | “Only six of 11 S. enterica isolates predicted to be azithromycin resistant were resistant.” |
| `38% of inconsistent calls (532/1,403)` | Two aminoglycosides dominated the errors — gene presence and phenotype correspond poorly for gentamicin and streptomycin in Salmonella | “accounting for 38% of inconsistent calls (532/1,403)” |
| `17%` | Per-isolate, rather than per-test, error rate: fraction of isolates with at least one genotype-phenotype mismatch | “In total, 17% of all isolates” |
| `21% of isolates tested against 12 antibiotics` | Why per-isolate errors are common even at high per-drug accuracy — the authors' own arithmetic | “21% of isolates tested against 12 antibiotics with a consistency of 98% would” |
| `< 90% identity` | Default BLASTP identity cutoff below which protein matches are discarded (per-gene curated cutoffs override it) | “(default of < 90%) or, if designated, a BlastRule cutoff are dropped” |
| `< 50% length of reference` | Default coverage cutoff below which partial matches are discarded | “(by default < 50% length of reference)” |
| `less than 90% length` | Definition of the PARTIAL/PARTIALP method label in the output | “The method PARTIAL is returned when a blast hit of less than 90% length is internal to a contig.” |
| `52 biocide, 8 heat, 148 metal` | Composition of the non-antibiotic (stress) part of the catalog that --plus turns on, database version 2020-07-16.2 | “resistance genes, 52 biocide resistance genes, 8 heat resistance genes, and 148 metal resistance genes.” |
| `34.2%` | Carrying no resistance genes is a common, normal result in the validation collection | “2,136 of the 6,242 isolates (34.2%) were pansusceptible” |

### Not established from the papers (5)

- Neither paper reports a sensitivity/specificity for the virulence or stress-gene calls
  themselves; the 2021 validation is two small Salmonella sets (19 and 6 isolates) compared
  against earlier surveys, not a systematic benchmark, so the --plus columns are far less
  validated than the AMR ones.

- No genotype-phenotype concordance is reported for any non-foodborne taxon, for anaerobes,
  or for MAGs — the four species tested (S. enterica, C. coli, C. jejuni, E. coli) are the
  entire evidence base.

- The database counts quoted here are from catalog version 2020-07-16.2; CompareM2 fetches
  the current catalog with `amrfinder -u`, and the papers give no way to know how much
  larger or differently composed that version is.

- I could not establish from the papers what AMRFinderPlus writes in the Class column for
  virulence genes that have no assigned class (only stx and intimin are described as
  receiving class/subclass), so it is unclear whether a genome with virulence hits will show
  them under a blank or 'NA' column in this report table, or not at all — the CompareM2 test
  genomes happened to contain no virulence hits.

- Neither paper reports how often a hit reported by protein-only mode would have been
  rejected or reclassified had nucleotide sequence been supplied as well, so the size of the
  accuracy cost of CompareM2's protein-only choice is unquantified.

---

## mlst
*Jolley KA, Bray JE, Maiden MCJ 2018, Wellcome Open Research 3:124, doi:10.12688/wellcomeopenres.14826.1*

MLST puts a short, globally shared name on each genome: a sequence type (ST), which is just
a label for the exact combination of alleles found at a handful of housekeeping genes. It
answers "has this isolate's genotype been seen and named before, and is it the same clone as
the one in that other lab's paper?" — not "how closely related are these two genomes?"

**How it works.** Each assembly is BLAST-searched against PubMLST's catalogue of known allele sequences for
the loci of a typing scheme; every locus is assigned the identifier of the allele it
matches, and that combination of allele numbers (the profile) is looked up in the scheme's
profile table. As the paper puts it for a genome queried against MLST, "individual allelic
matches will be identified along with the ST if the combination of alleles has been
previously defined" — so the ST exists only if a curator has already registered that
profile. The allele and profile identifiers are a curated nomenclature maintained in
PubMLST/BIGSdb, not a computed measurement.

### Reading the output

**Scheme.** Which PubMLST species/genus scheme was auto-detected for that genome. CompareM2 runs `mlst`
bare — no scheme is forced — so the detection is per genome and different rows may carry
different scheme names. PubMLST hosted databases for over 100 species or genera at the time
of the paper; a genome outside that coverage (many environmental MAGs, novel taxa) will show
'-' here and nothing useful in the other columns. Two rows are only comparable if this
column reads the same.

**Sequence type.** An index into that scheme's profile table, nothing more. It is scheme-relative: ST 5 in one
scheme has no relationship whatsoever to ST 5 in another, and the number carries no ordering
— ST 78 is not 'between' ST 77 and ST 79. A '-' here means the observed allele combination
was not in the profile table: either one or more loci failed to call, or the combination
itself is novel and would need to be submitted to PubMLST to be given a number.

**Alleles.** One entry per locus, e.g. `arcC(1), aroE(4), glpF(1), ...`. Conventional MLST schemes are
small — the paper describes MLST as indexing "multiple, but few (six or seven), housekeeping
gene fragments" — so expect roughly six or seven entries for a classical scheme. The number
in brackets is a catalogue identifier assigned when that sequence was first curated, not a
similarity score: allele 4 is not more similar to allele 5 than to allele 100.

**'-' or a marked allele.** A '-' at a locus means no allele was called there; a marked call (the tseemann/mlst software
flags inexact and partial matches with its own symbols next to the number) means the locus
matched a known allele only approximately or over part of its length, i.e. possibly a novel
allele. The paper is explicit that a genuinely core gene can be absent from a WGS dataset
because of "technical issues due to incomplete" assembly, so in a fragmented draft or MAG a
'-' is more often a contig break than real gene loss. A novel sequence needs curation before
it becomes an allele — BIGSdb allele definition "requires curator oversight if a new
sequence" is too different from an existing one. Any '-' or inexact locus normally means no
ST.

**Genome.** One row per input assembly, taken from a single `mlst` run over all of them; CompareM2 maps
the file path back to your sample name. Genomes that typed to no scheme still get a row, all
dashes — that is a real result (no applicable scheme), not a crash.

### Caveats

- An ST is meaningless without its scheme: the paper describes a scheme as an arbitrary
  grouping of loci chosen for a purpose, with no limit on how many schemes exist, so STs
  must never be compared across schemes, across species, or against STs from a different
  database — most schemes now live on PubMLST, but the paper notes that Salmonella and
  Escherichia coli are hosted on Enterobase and only mirrored there.

- A dash is usually about your assembly, not your organism: the paper states that core genes
  go missing from WGS datasets partly because of "technical issues due to incomplete"
  assembly, so fragmented drafts and MAGs lose loci at contig breaks and end up with no ST;
  check assembly contiguity before concluding a genome is untypeable or novel.

- Six or seven housekeeping genes is a coarse ruler — the paper's whole argument is that
  outbreak-level resolution needs cgMLST/rMLST/wgMLST (its cgMLST clustering thresholds go
  down to "5 or fewer" allele differences for point-source outbreaks, at schemes of
  ~1500-2000 loci). Two genomes sharing an ST are the same lineage, not the same isolate;
  use the ANI, SNP-distance or tree sections for relatedness.

- Not from the paper, but from how the tool is packaged: allele and ST numbers come from a
  static copy of PubMLST bundled with the installed `mlst` version, while PubMLST itself
  grows continuously with submissions — a profile with no ST today can acquire one after a
  database update, so record the tool/database version alongside any ST you publish.

### Quantitative record (7 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `six or seven` | Conventional MLST indexes only a small number of housekeeping gene fragments, so a classical scheme has roughly six or seven loci in the Alleles column | “multiple, but few (six or seven), housekeeping gene fragments” |
| `over 100` | Number of species/genera for which PubMLST hosted databases at the time of writing (2018) — the ceiling on which genomes can be typed at all | “databases for over 100 species or genera, mainly bacteria,” |
| `more than 300,000 isolate records and 100,000 genome assemblies` | Content of the PubMLST databases in 2018, i.e. the reference set your genome's profile is looked up against | “more than 300,000 submitted isolate records and 100,000” |
| `approximately 125 curators and 2000 active data submitters` | Allele and ST numbers are human-curated nomenclature, maintained by a curator community of this size | “approximately 125 curators and 2000 active data submitters” |
| `n=12,179 genomes, cgMLST with 1605 loci` | Seven-locus MLST clonal complexes agreed well with the finer cgMLST population structure in a large Neisseria meningitidis dataset — evidence that classical STs do capture lineages, even though they are coarse | “Neisseria meningitidis genomes (n=12,179) differentiated by cgMLST (1605 loci)” |
| `1998` | The MLST nomenclature your ST comes from has been accumulating since the first scheme, so ST numbering reflects order of discovery over decades | “the first MLST scheme developed in 1998” |
| `about a second` | Typing one assembly against a scheme is computationally trivial (measured for a query submitted to the PubMLST website, not for this local tool) | “The analysis usually takes about a second to perform.” |

### Not established from the papers (6)

- The paper gives no identity or coverage cut-off for calling a seven-locus MLST allele from
  a draft assembly. Its numeric rules (new cgMLST/pan-genome alleles must be "within 98%
  identity and 98% total length of a known allele"; exemplar alleles within 10% identity)
  describe BIGSdb's own curation and BLAST search-space heuristics, not what tseemann/mlst
  does. Whatever thresholds that tool applies are undocumented in this paper, and CompareM2
  passes no arguments that would change them.

- The exact output legend — which characters mark an inexact allele, a partial allele, a
  missing allele and multiple hits at one locus — is defined by tseemann/mlst, which has no
  publication. `mlst` is not installed in this environment, so I could not confirm the
  legend against a real run; the report currently prints those symbols through verbatim with
  no key.

- No threshold is given anywhere in the paper for how many differing loci in a seven-locus
  profile constitute 'related' or the same clonal complex; the clustering thresholds it does
  give (200 or more for major lineages, 5 or fewer for point-source outbreaks) are for
  cgMLST schemes with hundreds to thousands of loci.

- Scheme auto-detection is not evaluated in the paper, and the output carries no confidence
  value for it. If two genomes in your set carry different scheme names, nothing in the
  table tells you whether that is real biology or a detection artefact.

- The paper mentions tseemann/mlst only as a table entry among third-party tools using
  PubMLST data; it neither describes nor benchmarks its algorithm, so nothing here is a
  validation of the specific implementation CompareM2 runs.

- CompareM2 v3 runs `mlst <assemblies>` with no arguments, so there is no way to force a
  scheme from the config in this version (v2 exposed `set_mlst--scheme`) and the list of
  available schemes is no longer written out.

---

## mashtree
*Katz LS, Griswold T, Morrison SS, Caravas JA, Zhang S, den Bakker HC, Deng X, Carleton HA 2019, Journal of Open Source Software 4(44):1762, doi:10.21105/joss.01762 — underlying method: Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM 2016, Genome Biology 17:132, doi:10.1186/s13059-016-0997-x*

Mashtree draws a tree of all your genomes at once without aligning anything. It compresses
each assembly into a small fingerprint of hashed 21-mers, turns the fraction of fingerprints
two genomes share into a distance that tracks 1 minus their average nucleotide identity, and
joins everything into a tree from that distance matrix — so it answers "which of my genomes
group together, roughly", in seconds rather than hours.

**How it works.** Mash slides a window of length k=21 over each assembly, hashes every canonical k-mer, and
keeps only the 10,000 smallest hash values as that genome's "bottom sketch"; the fraction of
hashes two sketches share is an unbiased estimate of the Jaccard index, which Mash converts
into the Mash distance D = -(1/k) ln(2j/(1+j)), an estimate of the point-substitution rate
under a simple Poisson model. Mashtree feeds the resulting all-against-all Mash distance
matrix to QuickTree's neighbour-joining implementation with default options to produce a
dendrogram, which CompareM2 writes to `<output_directory>/mashtree/mashtree.newick` and
renders as the phylogram in this section.

### Reading the output

**Horizontal branch lengths in the tree.** Distance along the horizontal axis is Mash distance, which the Mash paper shows tracks 1 -
ANI (D ~= 1 - ANI) with a root-mean-square error of 0.00274 against alignment-based ANI at
k=21 for 500 Escherichia genomes. As a rule of thumb from the same paper, two genomes
separated by a total path of about 0.05 sit at roughly 95% ANI, the level used there to
build species-level clusters. Vertical spacing is only there to keep leaf labels apart and
carries no information.

**No scale bar and no axis.** The report draws the phylogram scaled to fit the page, with no ruler. You can compare branch
lengths to each other, but you cannot read an absolute Mash distance off the picture — for
that, open `<output_directory>/mashtree/mashtree.newick`, where the branch lengths are
printed as numbers.

**Absence of node support values.** There are no percentages on the internal nodes, and this is not a rendering omission.
Mashtree does implement bootstrapping and jackknifing, but CompareM2 runs plain `mashtree`,
so no resampling is done and the newick carries no support values. Treat every split as
unmeasured: a group that looks tight here has not been tested for stability.

**The leftmost node (apparent root).** Neighbour-joining produces an unrooted tree, so the node the drawing starts from on the left
is a drawing convention, not an inferred common ancestor. Do not read "earliest-branching"
or "most ancestral" off the left edge — only the grouping and the distances between leaves
are meaningful. (That NJ output is unrooted is standard method knowledge; the Mashtree paper
simply calls the output a dendrogram and states that Mashtree does not infer phylogeny.)

**Leaf labels with a " · <name>" suffix and colour.** When TreeCluster has also run, each genome name is followed by its cluster assignment and
coloured by it, with the note "Leaves coloured by TreeCluster assignment." above the figure.
Those clusters are cut from this same tree by a separate tool at a branch-length threshold;
they are a summary of this drawing, not independent evidence about it.

### Caveats

- This is not a phylogeny and the authors say so: Katz et al. write "Although Mashtree does
  not infer phylogeny", and Ondov et al. "emphasize that Mash is not explicitly designed for
  phylogeny reconstruction, especially for genomes with high divergence or large size
  differences" — use it to triage and group, and use CompareM2's core-genome tree
  (fasttree/IQ-TREE on the Panaroo alignment) for any evolutionary claim.

- CompareM2 runs plain `mashtree`, not its bootstrap or jackknife workflows, so no node in
  this tree has a support value; the confidence-value machinery described in the Mashtree
  paper is simply not exercised here.

- The Mash distance measures resemblance of whole k-mer sets, so it moves with gene content
  and genome size as well as with point mutations: the Mash paper notes that J is sensitive
  to genome size and that Mash deliberately sets the denominator to the average genome size
  to penalise size differences. A draft genome that is more fragmented, more contaminated,
  or carries extra plasmids will therefore shift on this tree for reasons that are not
  divergence — cross-check against the CheckM2 completeness/contamination and seqkit
  assembly-size columns before believing an odd placement. (The connection to assembly
  quality is general caution; the papers state the genome-size sensitivity, not this
  specific failure mode.)

- The demonstrated accuracy is for closely related genomes — strong correlation with ANI was
  shown over 90–100 % ANI, and the Mash paper states the correlation degrades for more
  divergent genomes because the variance of the estimate grows with distance. Deep splits
  near the base of a mixed-genus set are the least trustworthy part of the drawing.

### Quantitative record (13 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `1,000 to 10,000` | Mashtree does not use Mash's default sketch size; it raises the number of hashed k-mers kept per genome tenfold. CompareM2 passes --sketch-size 10000 explicitly, matching Mashtree's own choice. | “from 1,000 to 10,000 to increase discriminatory power” |
| `D ≈ 1 − ANI` | The Mash distance is designed to approximate 1 - ANI, and does so across several sketch and k-mer settings. | “D ≈ 1 − ANI over multiple sketch and k-mer sizes” |
| `0.00274` | Accuracy of the Mash distance as a stand-in for 1 - ANI, measured against MUMmer dnadiff ANI on 500 Escherichia genomes at k=21 with the smaller default sketch of 1000. | “a root-mean-square error of 0.00274” |
| `90–100 %` | Range of divergence over which the ANI/Mash-distance correlation was demonstrated to be strong. | “For ANI in the range of 90–100 %, the correlation with Mash” |
| `≤0.05` | The Mash-distance cut-off the authors used to form species-level connected components across all of RefSeq. | “RefSeq genomes with a pairwise Mash distance ≤0.05,” |
| `≥95 %` | The ANI level that Mash-distance cut-off corresponds to. | “which equates to an ANI of ≥95 %” |
| `k = 14 and k = 19` | k=21 (CompareM2's --kmerlength) is well above the k needed to avoid chance k-mer collisions at CompareM2's declared --genomesize of 5,000,000, and also above what a 3 Gbp genome would need, so spurious k-mer sharing is not a concern for bacterial or archaeal assemblies. | “which yields k = 14 and k = 19 for 5 Mbp and 3 Gbp genomes (q = 0.01), respectively” |
| `k = 21` | k=21 is the setting the Mash authors found accurate in most cases and set as the default; CompareM2 keeps it. | “We have found the parameters k = 21 and s = 1000 give accurate estimates in” |
| `0.10` | How far a Mash-distance neighbour-joining tree drifted from an alignment-based tree in the only tree benchmark either paper reports (17 primate genomes vs the UCSC alignment tree), with Mash branch lengths slightly longer on average. | “Branch Score Distance [29] of 0.10 between the two” |
| `1000 to 5000` | Too small a sketch can misplace a taxon outright: adding five more divergent mammals to the primate set misplaced the tarsier at sketch size 1000, and enlarging the sketch fixed it. This is the direct evidence that sketch size, not just k, controls topology at larger distances. | “Increasing the sketch size from 1000 to 5000 corrects” |
| `ANI >95 %` | For genomes as similar as those in a typical CompareM2 run, even a few hundred hashes suffice for basic clustering — CompareM2's 10,000 is comfortably above that; the sketch size matters mainly for divergent pairs. | “(e.g. ANI >95 %), sketches of a few hundred hashes are” |
| `0.99` | Confidence level of the published Mash-distance error bounds (Table 1), which shrink as sketch size grows. | “For a given sketch size and Mash distance, the Mash estimation error will be less than the given value with 0.99 probability” |
| `54,118 genomes in 33 CPU h` | The scale at which this approach is intended to work, for context on why it is in the pipeline at all. | “the clustering of all 54,118 NCBI RefSeq genomes in 33 CPU h” |

### Not established from the papers (6)

- Neither paper documents what mashtree's --mindepth 5 (which CompareM2 passes) does for
  assembled FASTA input. The Mash paper only describes low-abundance k-mer filtering in the
  context of raw sequencing reads, where it assumes redundancy at "depth of coverage >5",
  and notes that by default the coverage threshold is set to one; I could not establish from
  the papers whether this flag has any effect on the assemblies CompareM2 supplies.

- Likewise, neither paper says what --genomesize 5000000 does for assembled input; in Mash
  the genome-size hint is discussed in relation to Bloom-filter memory when sketching reads
  and in relation to choosing k, not as something that changes the distance for a finished
  assembly.

- No Mash-distance or branch-length threshold for anything below the species level is given
  in either paper. The only threshold stated is D ≤0.05 / ANI ≥95 % for species-level
  connected components, and even that is tied to the 70 % DNA-DNA reassociation definition
  that the paper itself calls "a historical, albeit debatable, definition of bacterial
  species". CompareM2's TreeCluster default of 0.05 is a pipeline choice, not a value either
  paper validates for max-clade tree cutting.

- Neither paper benchmarks mashtree topology against a reference bacterial or archaeal
  phylogeny. The only tree comparison is 17 primate genomes against a UCSC alignment tree —
  eukaryotic, and with n=17.

- The Mashtree paper reports no runtime, memory or accuracy benchmarks at all, so I cannot
  give an expected wall-clock cost for this step in CompareM2.

- The papers do not state whether QuickTree's newick output is emitted rooted or unrooted,
  nor how negative branch lengths (a known NJ artefact) are handled; CompareM2's renderer
  would draw a negative length as a backwards branch.

---

## treecluster
*Balaban M, Moshiri N, Mai U, Jia X, Mirarab S 2019, PLoS ONE 14(8):e0221068, doi:10.1371/journal.pone.0221068*

TreeCluster takes the mashtree that CompareM2 already built and chops it into groups, so
each genome gets a cluster number instead of you having to eyeball the tree. Two genomes end
up in the same cluster only if they sit close together on that tree — here, within 0.05 of
mash distance.

**How it works.** TreeCluster solves a "min-cut partitioning" problem: given a tree and a threshold, cut the
fewest edges such that every resulting group of leaves stays under a diversity limit.
CompareM2 runs `TreeCluster.py -i mashtree.newick --method max_clade --threshold 0.05`,
which sets the limit to the longest leaf-to-leaf path inside a group (the group's diameter)
and additionally requires each group to be a clade — a node of the tree plus all of its
descendants — which the paper solves by cutting both child edges at a node instead of only
the longer one.

### Reading the output

**SequenceName (first column of treecluster.tsv).** The leaf label copied straight from mashtree.newick, i.e. one row per input assembly. The
treecluster section of the report is a plain two-column dump of this file, so the names are
whatever mashtree wrote, not a tidied-up sample list.

**ClusterNumber (second column).** An arbitrary integer label — 1, 2, 3 ... — with no ordering or meaning beyond "these genomes
were grouped together". Each numbered cluster is a clade of the mashtree whose widest
internal leaf-to-leaf distance is at most 0.05. That last part is what --method max_clade
buys you over the alternatives: plain `max` enforces the same diameter limit but lets a
cluster be any connected piece of the tree rather than a whole clade; `sum` caps the total
branch length inside a cluster instead of its widest pair; `single_linkage` only requires a
chain of short hops, which the paper shows lets a very heterogeneous set collapse into one
cluster through transitivity.

**ClusterNumber = -1.** Not an error and not a failure code: TreeCluster writes -1 for every genome that ended up
alone, so all -1 rows are separate singletons, not one big cluster. Singletons are routine —
on the paper's 16S benchmark at threshold 0.09, 27% of TreeCluster's clusters were
singletons.

**"Leaves coloured by TreeCluster assignment" (mashtree section).** When treecluster.tsv exists, the mashtree phylogram colours each leaf by its cluster and
appends " · <cluster>" to the label. A leaf that stays in the default colour with no suffix
means its name in the tree did not match any row in treecluster.tsv, which is a naming
mismatch rather than a biological result.

**"First 50 rows.".** The generic table renderer truncates at 50 rows, so with more than 50 genomes the report
shows a partial assignment list; the complete table is treecluster/treecluster.tsv in the
output directory.

### Caveats

- The 0.05 threshold is in the branch-length units of the tree it was handed — mash
  distances from mashtree — so it means "no two genomes in a cluster are more than 0.05
  apart along the tree", not 95% ANI, not 97% identity, and not a species boundary; the
  paper is explicit that similarity thresholds are used in the first place to avoid defining
  species ("The use of a similarity threshold instead of a biological concept of species is
  to avoid the notoriously difficult problem of defining species for microbial organisms").

- Clusters are an artefact of the threshold, not a discovered structure: the paper's own
  sweep moved the answer from 181,574 to 10,112 clusters, so before reporting a grouping,
  re-run TreeCluster.py at a couple of neighbouring thresholds and see which genome pairs
  stay together.

- max_clade requires every cluster to be a clade of the tree as rooted in the newick file,
  and the paper's clade-constrained experiments rooted their trees first with
  MinVar/FastRoot; CompareM2 passes mashtree's newick through unmodified, so the rooting is
  whatever mashtree wrote — the paper does not measure how much an arbitrary rooting changes
  clade-constrained results, so treat this as general caution rather than a demonstrated
  problem.

- The optimal clustering is not unique — the paper shows the number of equally optimal
  partitions can be exponential in the number of leaves — so for a genome sitting right at
  the threshold boundary, its specific cluster assignment is one of several equally valid
  answers.

- TreeCluster never looks at the genomes; it only reads the tree, so any error in the mash
  sketching or the neighbour-joining topology is inherited whole, and the clusters cannot be
  more reliable than the mashtree they were cut from.

### Quantitative record (7 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `181,574 and 10,112 clusters` | Cluster membership is strongly threshold-dependent: on the Greengenes tree of 203,452 sequences, sweeping the threshold from 0.005 to 0.15 changed the number of clusters by roughly 18-fold. | “181, 574 and 10, 112 clusters (note that singletons are also counted).” |
| `20 thresholds from 0.005 to 0.15` | The thresholds the authors themselves explored for OTU clustering on tree distances span a wide range; there is no single recommended value. | “use the following 20 thresholds: [0.005, 0.05] with a step size of 0.005, and (0.05, 0.15] with a” |
| `4.5%` | The threshold the paper used for Max-diameter clustering in its HIV application was inherited from another tool's default, i.e. it is application-specific and not a property of the method. | “threshold for Max-diameter clustering is 4.5% [4], so we use this as our clustering threshold” |
| `27%` | A large fraction of TreeCluster output rows can be singletons — in the paper's 16S benchmark at threshold 0.09, just over a quarter of clusters contained a single sequence. | “27% of the clusters are singletons for TreeCluster.” |
| `22,090 and 23,631` | At threshold 0.09 on the Greengenes tree, TreeCluster's Max-diameter mode produced 23,631 clusters, comparable in number to the UCLUST-based Greengenes OTUs. | “for α = 0.09, both methods have similar number of clusters (22,090 and 23,631 for” |
| `under 2 seconds for 5,000 leaves` | The clustering step itself is negligible in runtime once the tree exists, so re-running it at several thresholds is cheap. | “set with 5,000 leaves, the running time of TreeCluster did not exceed 2 seconds.” |
| `30 seconds for more than 200,000 leaves` | Even on a very large tree the clustering finishes in well under a minute. | “leaves, TreeCluster performed clustering in only 30 seconds.” |

### Not established from the papers (5)

- The paper gives no recommended threshold for whole-genome mash-distance trees. Its
  thresholds (0.005-0.15 on a 16S ML tree, 4.5% on HIV pol) come from different data types
  on different distance scales; CompareM2's 0.05 is a pipeline default and I could not trace
  it to anything in this paper.

- TreeCluster is never evaluated on whole-genome or mash-distance data anywhere in the paper
  — the three applications are 16S OTU clustering, HIV transmission clustering, and
  multiple-sequence-alignment decomposition. Its performance on a mashtree of bacterial
  genomes is untested here.

- The paper compares max_clade to unconstrained max only indirectly: it reports that in the
  HIV simulation the "Clade constraint has little impact on effectiveness", and gives no
  such comparison for OTU clustering. Whether the clade constraint changes anything on a
  genome tree is not established.

- I could not determine what the number of clusters or their sizes should look like for a
  well-behaved set of genomes; the paper offers no diagnostic for a bad clustering other
  than sweeping the threshold.

- Whether mashtree's newick is midpoint-rooted or written with an arbitrary root is a
  property of mashtree, not of this paper, and I did not verify it.

---

## skani
*Shaw J & Yu YW 2023, Nature Methods 20:1661–1665, doi:10.1038/s41592-023-02018-3*

skani answers "how similar is each of my genomes to each of the others?" by computing
average nucleotide identity (ANI) — the percent identity over the parts of two genomes that
share ancestry. The report shows it as a shaded genome-by-genome matrix, which is the
fastest way to spot duplicates, see which inputs are the same species, and find the odd one
out.

**How it works.** skani takes a very sparse subset of k-mers from each genome (defaults k = 15, roughly one
k-mer in c = 125), chains matching k-mers to locate orthologous regions, cuts the query
genome into 20-kb chunks, and estimates identity per chunk from the fraction of its seeds
that anchor into a chain, then takes a weighted mean over chunks and applies a learned
regression that debiases the estimate against a MUMmer-based (ANIm) reference. Because
identity is measured only inside regions that actually chain, sequence missing from an
incomplete assembly does not drag the ANI down the way it does for pure sketching methods
such as Mash.

### Reading the output

**Matrix cells (ANI, %).** Each cell is the estimated percent identity over the regions of the two genomes that could
be matched. The diagonal is 100.00 (a genome against itself), and the matrix is symmetric:
skani chooses which genome is 'reference' by size, so the value "does not depend on the
order of the inputs" and the two halves of the full matrix mirror each other. The paper
refers to "the standard 95% ANI species threshold" (citing Jain et al. 2018) as the
convention for calling two genomes the same species; skani's estimator is described as most
accurate at ANI >= ~82%, so this matrix is a within-species and near-species instrument.

**"Shaded from X% to 100%" (the line above the table).** The colour scale is stretched to the lowest value actually present in your run, not to a
fixed scale. If all your genomes are one species, the scale may run 98-100% and tiny
differences will look dramatic; if one distant genome is in the set, the scale bottom drops
and everything else washes out to the same shade. Read the numbers; the colour only ranks
them within this particular run.

**Empty or near-zero cells.** skani refuses to report a pair in two situations: the fast marker screen puts putative ANI
below 80% (only pairs with ANIFMH > 80% are compared), or the aligned fraction is under
~15%. A cell with no real value therefore means 'too little detectable homology to
estimate', not '0% identical'. Pairs from different genera normally land here. If such a
cell is written as a literal 0.00 it will also drag the shading scale down to zero, so check
the printed range line if the whole table suddenly looks pale.

**Aligned fraction (AF) — computed by skani but not shown in this table.** AF is the share of each genome covered by the matched regions (the paper defines it as "the
sum of the bases covered by" all chains divided by genome size). CompareM2 runs `skani
triangle --full-matrix`, which writes an ANI matrix only, so the report shows identity
without telling you how much of each genome that identity covers. Since an ANI is emitted
with AF as low as ~15%, a high number between a full chromosome and a small plasmid or a
partial MAG is possible and means something quite different from high identity across both
whole genomes. Check the assembly sizes in the QC section before reading a high ANI as 'same
organism'.

**Differences within the 98-100% band.** Against the BLAST-style OrthoANIu baseline on 4,350 E. coli genomes, skani had Pearson R and
mean absolute error of (0.981, 0.131) — so gaps of one or two tenths of an ANI point are
inside the method's own error and should not be read as real. On the diverse all-to-all
refseq-rc set, which spans much lower identities, the same statistics were (0.976, 1.279):
precision degrades substantially once pairs cross species boundaries.

### Caveats

- ANI without AF is half the picture: skani emits a value once the aligned fraction reaches
  only about 15% of one genome, and the matrix in this report carries no AF column, so a
  high identity between genomes of very different size or completeness can reflect a shared
  plasmid, prophage or conserved core rather than whole-genome relatedness.

- This matrix cannot resolve anything below species level: skani screens out pairs whose
  putative ANI is under 80%, is accurate only at ANI >= ~82%, and the authors list low ANI
  among the extreme regimes where skani is limited — so an absent or zero cell means 'more
  distant than skani measures', and gives no information about how distant.

- CompareM2 runs skani at defaults (k = 15, c = 125), which the paper describes as tuned for
  similar genomes; for fragmented assemblies the authors recommend c = 70 (ANI <= 95, N50 <=
  10 kb) or c = 30 (N50 ~ 3 kb) and report that lowering c improves both ANI and aligned-
  fraction accuracy, so MAG-heavy genome sets are being run at the least accurate of the
  described settings.

- The 95% figure is a convention the paper cites, not one it validates, and skani's own
  values are only regression-corrected toward an alignment-based reference above 90%
  putative ANI; a pair landing at 94.8 or 95.2 should not settle a species assignment on
  this table alone (general methodological caution, beyond what the paper claims).

### Quantitative record (18 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `ANI >= ~82%` | ANI range in which skani's estimator is accurate | “most accurate when ANI ≥ ~82%” |
| `> 80%` | Sketch-based pre-screen: pairs below this putative ANI are never compared, so no value is emitted for them | “ANIFMH > 80% as a conservative underestimate.” |
| `AF >= 15%` | Minimum predicted aligned fraction for skani to output an ANI value at all (default) | “predicted AF ≥15% by default, which ends up giving reasonable ANIs” |
| `82%` | Lowest ANI at which the authors consider the values reasonable on their three benchmark datasets | “down to the 82% range on the three datasets shown.” |
| `95% ANI` | The species-delineation convention the paper refers to (citing Jain et al. 2018); the paper cites it, it does not establish it | “when subject to the standard 95% ANI species threshold” |
| `R2 = 0.001` | Dependence of skani's ANI on assembly incompleteness and contamination among MAG pairs that are >99% identical by ANIm (ordinary least-squares R2; near zero means the estimate is not distorted by assembly quality) | “R2 = 0.001” |
| `41,494 comparisons` | Number of MAG pairs in that incompleteness/contamination benchmark | “Skani (41,494 comparisons)” |
| `4% ANI at 50% completeness` | How far a pure sketching method (Mash) can be pulled down by incompleteness alone — the reason skani measures identity only inside aligned regions | “to a difference of 4% ANI at 50% completeness” |
| `(0.981, 0.131)` | skani versus the BLAST-style OrthoANIu baseline on a collection of E. coli genomes (Pearson R, mean absolute error in ANI points) | “(0.981, 0.131)” |
| `n = 4,350` | Size of that E. coli benchmark | “n = 4,350” |
| `(0.976, 1.279)` | Same comparison against OrthoANIu on the diverse refseq-rc all-to-all set, which spans lower identities (Pearson R, mean absolute error in ANI points) | “(0.976, 1.279)” |
| `n = 4,233, m = 4,233` | Size of the refseq-rc all-to-all benchmark | “n = 4,233, m = 4,233” |
| `>500 times faster than ANIm and >50 times faster than FastANI` | Speed of computing an all-against-all distance matrix relative to alignment-based ANIm and to FastANI | “is >500 times faster than ANIm and >50 times faster than FastANI for” |
| `putative ANI > 90%` | Range over which the learned debiasing correction is applied; below it the printed values are raw chaining estimates | “We only debias comparisons with putative ANI > 90%” |
| `k = 15 and c = 125` | Default seeding parameters, which is what CompareM2 runs (k-mer length; one k-mer in c is kept as a seed) | “By default, k = 15 and c = 125.” |
| `20-kb nonoverlapping chunks` | Chunk size the query genome is cut into before per-chunk identity estimation | “fragment the query into 20-kb nonoverlapping chunks” |
| `c = 70` | Non-default setting the authors recommend for genomes below 95% ANI or with N50 at or below 10 kb | “with c = 70 for genomes with ANI ≤95 and N50 ≤10 kb” |
| `c = 30` | Setting the authors recommend for the most accurate aligned-fraction estimates on very fragmented genomes | “empirical heuristics: a ‘slow’ pre-set with c = 30 for the most accurate” |

### Not established from the papers (4)

- What skani literally writes into a matrix cell for a pair it filtered out (a blank, or
  0.00) is not stated in the paper, and skani is not installed in this environment, so I
  could not check; the report prints whatever is in ani.tsv and would also stretch its
  colour scale down to that value.

- The paper gives no guidance on how much aligned fraction is 'enough' for a comparison to
  be trusted beyond the 15% output gate, and no AF-based threshold for calling two genomes
  the same organism.

- No ANI bands are given for genus, family or any rank below species, so the matrix cannot
  be annotated with taxonomy beyond the cited 95% species convention.

- Behaviour on plasmid-only, phage or very small contig sets is not addressed; the non-
  prokaryotic benchmarking is on eukaryotic MAGs with a median size of 17.6 Mbp, i.e. the
  opposite end of the size range.

---

## panaroo
*Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. 2020. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biology 21:180. doi:10.1186/s13059-020-02090-4. Supporting components: Fu L, Niu B, Zhu Z, Wu S, Li W. 2012. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics 28(23):3150-3152, doi:10.1093/bioinformatics/bts565; Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7. Molecular Biology and Evolution 30(4):772-780, doi:10.1093/molbev/mst010.*

Panaroo sorts every predicted gene in your genomes into gene clusters (roughly: gene
families) and tells you which clusters are in all of your genomes, which are in some, and
which are in only one. It answers "what gene content do these genomes share, and what makes
each one different?" — and it does so while actively trying to undo the annotation mistakes
that otherwise make that answer wrong.

**How it works.** Panaroo clusters all annotated gene sequences with CD-HIT at 98% identity, then builds a
graph in which each node is a cluster of orthologous genes and an edge joins two nodes if
those genes sit next to each other on a contig in any genome; because a genome may appear
only once per cluster, paralogues are split out and then re-collapsed using that
neighbourhood context. It then uses the graph to clean up annotation error — merging genes
that are really one sequence translated in different frames, collapsing diverse families
that were over-split (re-compared at 70% identity when they share a neighbour), recursively
deleting poorly supported degree-1 nodes (contig-end artefacts and contamination), and re-
searching the DNA around a neighbouring node for genes the annotator simply missed.
CompareM2 runs `--clean-mode strict` (the aggressive deletion setting) with `-a core`, which
additionally writes a multiple-sequence alignment of the core clusters that snp-dists and
FastTree then consume.

### Reading the output

**Summary line: "N gene clusters across M genomes, in K distinct presence patterns".** N is the pan genome — every gene cluster seen anywhere in this set. K is how many different
on/off combinations across your genomes actually occur; the report notes the hard ceiling is
2**M - 1. K close to 1 means the genomes are near-identical in gene content; K close to its
ceiling means gene presence is scattered rather than tracking a few lineages. Judge N
against expectation for the taxon: in the paper's control dataset of 413 clonal M.
tuberculosis genomes the true answer was a core of roughly 4000 genes and essentially no
accessory, and tools that got it wrong got it wrong by 2584-3670 spurious accessory genes.

**The presence matrix figure (core on the left, genome-specific on the right).** One row per genome; each vertical block is one presence pattern, its width proportional to
how many clusters share that pattern. The solid block at the far left, filled in every row,
is the core. Ragged blocks in the middle are shared by subsets — these often follow
phylogeny, or a shared plasmid or prophage. Blocks at the far right touch one row only:
genes unique to that genome. A single row that is visibly emptier than the others across the
whole figure usually means a fragmented or incomplete assembly rather than real gene loss.

**"Shared gene content" table — Pattern / Present in / Gene clusters / Share.** Ranked by how many clusters carry each pattern. The top row is almost always "all genomes"
(the core) with the largest Share. Rows labelled "only <genome>" are that genome's private
gene content; under strict clean-mode these have already survived aggressive pruning, and
the paper's rationale is that this "retains rare genes that have reliable contextual
support" while deleting the rest. A row shared by an unexpected subset of genomes is the
interesting one — check whether those genomes are also neighbours in the tree; if they are
not, suspect horizontal transfer or a shared mobile element.

**"Pangenome partitions" — Core (>=99%) / Soft core (95-99%) / Shell (15-95%) / Cloud (<15%).** The 99% core cutoff is the convention the Panaroo paper uses, described there as "the 99%
presence threshold for core genes as used in Roary". Note the arithmetic of these bins in a
small run: with fewer than 100 genomes, >=99% can only be reached by a cluster present in
every genome, and a cluster missing from exactly one of M genomes sits at (M-1)/M, which
only reaches 0.95 once M >= 20. So with a typical CompareM2 set of a handful of genomes,
Soft core is 0 by construction and one missing gene call drops a cluster straight from Core
to Shell. Treat a slightly small Core with a fat Shell as a possible annotation artefact,
not automatically as biology.

**"Genes per genome" — Genome / Gene clusters / Of pangenome.** Gene clusters is that genome's total (core plus whatever accessory it carries); "Of
pangenome" is its share of the N clusters in the whole set. Genomes of the same species
should sit within a few percent of each other. One genome markedly below the rest points at
an incomplete assembly or missed annotations; one markedly above points at contamination or
a genuinely larger accessory complement (plasmids). Cross-check against the CheckM2
completeness and contamination values before concluding anything biological — the paper
recommends this kind of pre-screening explicitly.

### Caveats

- CompareM2 runs `--clean-mode strict`, which recursively deletes poorly supported degree-1
  nodes, and the paper states plainly that this "can occasionally lead to rare plasmids
  being removed" — so an absent gene in this matrix is not evidence of absence, and this run
  cannot be used to argue that a strain lacks a rare plasmid.

- The paper says Panaroo "is not recommended for metagenomic datasets", so a pangenome built
  from MAGs inherits both that warning and the MAGs' own incompleteness, which will show up
  as an artificially small core rather than as real gene loss.

- Core and accessory here are defined only against the genomes in this run — the paper's own
  definition of a pangenome is "the set of all genes that have been found in a species as a
  whole", so adding or dropping one genome rewrites every partition count, and a handful of
  isolates does not estimate a species pangenome.

- The paper validated Panaroo on Prokka GFF3 input and recommends running its own QC script
  on the assemblies beforehand; CompareM2 feeds Bakta annotations and applies no such pre-
  filter, so a contaminated or fragmented assembly will still be pulled into the graph —
  read the CheckM2 section before trusting an unusual pan genome size (not a claim the paper
  makes about Bakta specifically).

### Quantitative record (11 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `98%` | Initial clustering of all gene sequences is done by CD-HIT at a high nucleotide/protein identity threshold | “high sequence identity threshold (98%)” |
| `70%` | Clusters that share a neighbour in the graph are re-compared at a relaxed threshold so that over-split diverse gene families get merged back together | “compared at a lower pairwise sequence threshold (default 70%)” |
| `95% coverage and 99% identity` | Coverage and identity at which two nearby genes are called a mistranslation/frameshift of the same gene and collapsed into one node | “typically 95% and 99% respectively” |
| `maximum pairwise SNP distance of 9` | In the M. tuberculosis control dataset the isolates were so clonal that essentially no accessory genome should exist | “the maximum pairwise SNP distance within this dataset was 9” |
| `2584 to 3670 genes, a nearly tenfold increase over Panaroo` | On that clonal dataset, the other pangenome tools reported accessory genomes that should not have existed at all | “2584 to 3670 genes representing a nearly tenfold increase to that reported by Panaroo” |
| `59%` | Share of the between-tool difference on the M. tuberculosis dataset attributable to genes being broken across contigs | “genes being fragmented during assembly (59%” |
| `12% of the dataset lost` | Cost of the alternative strategy — filtering out assemblies flagged by CheckM instead of correcting annotations — on that dataset | “in a loss of 12% of the dataset which could potentially have a large impact” |
| `3372 and 3376 genes` | Core genome size Panaroo inferred for 328 globally sourced Klebsiella pneumoniae genomes, in default and sensitive mode | “sensitive modes, 3372 and 3376 respectively” |
| `1800 genes` | Core genome size Roary inferred on the same K. pneumoniae dataset — the smallest of any method, attributed to clusters being incorrectly split | “Roary identified the smallest core genome of 1800 genes” |
| `99%` | Prevalence cutoff used in the paper to call a cluster 'core' — the same cutoff the CompareM2 report's Core row uses | “using the 99% presence threshold for core genes as used in Roary” |
| `95%` | Fixed identity threshold that made an older tool over-split gene families on a diverse dataset — the failure mode Panaroo's context-based 70% re-comparison is designed to avoid | “default Roary pairwise identity threshold of 95% is too stringent” |

### Not established from the papers (6)

- The prevalence cutoff Panaroo uses internally to decide which clusters enter
  core_gene_alignment.aln under `-a core` is not stated in the paper; the 99% figure quoted
  above is described there as the threshold used for reporting core gene counts, and I could
  not confirm the two are the same.

- The Soft core (95-99%), Shell (15-95%) and Cloud (<15%) boundaries in the report's
  partition table come from the report, not from the Panaroo paper, which gives only the 99%
  core threshold and no values for the other bins.

- The paper does not give a numeric node-support threshold for strict mode — it says only
  that Panaroo "recursively removes nodes of degree 1 that are" below a given support
  threshold — so I cannot quantify how much more aggressive strict is than sensitive, only
  that sensitive deletes nothing.

- The paper lists MAFFT, Prank and Clustal Omega as the alignment options and does not say
  which is the default; CompareM2 does not set `--aligner`, so the aligner actually used was
  not established from the papers.

- The report renders no pan-genome accumulation curve, which is consistent with the paper's
  warning that such curves "are not robust to errors and fail to account for sampling biases
  and population structure" and that "accumulation curves should not be used to compare
  pangenome characteristics of different lineages or species" — but that means the report
  gives no read on whether the pangenome is open or closed, and the paper's recommended
  alternatives (the IMG/FMG models, which need a dated phylogeny) are not run here.

- No runtime or memory expectation for a CompareM2-sized run: the paper's benchmark is for
  10, 100 and 1000 N. gonorrhoeae isolates with 5 CPUs and the values sit in a figure I did
  not read numerically.

---

## snp-dists
*Core gene alignment: Tonkin-Hill G, MacAlasdair N, Ruis C, et al. "Producing polished prokaryotic pangenomes with the Panaroo pipeline." Genome Biology 2020, 21:180, doi:10.1186/s13059-020-02090-4. snp-dists itself has no publication — it is cited in the CompareM2 paper only as software: Seemann T, Klötzl F, Page AJ. "Pairwise SNP distance matrix from a FASTA sequence alignment." tseemann/snp-dists GitHub, 2025. Pipeline: Kobel CM, Aho VTE, Øyås O, et al. Bioinformatics 2025, 41(9):btaf517, doi:10.1093/bioinformatics/btaf517.*

snp-dists counts, for every pair of your genomes, how many positions differ in the core gene
alignment that Panaroo built — the genes shared by (nearly) all genomes in this run. The
result is a square table of whole numbers: 0 means the two genomes are identical across that
shared core, and bigger numbers mean more single-base differences.

**How it works.** Panaroo clusters the annotated genes into orthologue groups, aligns the core ones and
concatenates them into a single alignment ("core and accessory genome alignments created
using either MAFFT, Prank or Clustal Omega"); snp-dists then walks that alignment column by
column and, for each pair of sequences, counts the columns where the two bases differ.
CompareM2 runs it as `snp-dists results/panaroo/core_gene_alignment.aln` with no options at
all, so what you get is the plain default: a symmetric matrix of raw, uncorrected difference
counts — no evolutionary model, no distance transformation, no normalisation by alignment
length.

### Reading the output

**Top-left header cell ("snp-dists 1.2.0").** This is not a column of data — snp-dists puts its own version string in the first header
cell, and the report shows the file as-is. The remaining column headers, and the first cell
of every row, are your genome names.

**The diagonal.** Every genome against itself is 0. A 0 off the diagonal means the two genomes are
indistinguishable across the core alignment, which is what you see for a genome included
twice, and also for genuinely identical isolates.

**Off-diagonal integers.** A raw count of differing alignment columns over the whole concatenated core alignment, not a
percentage and not a per-Mb rate. To compare a value against anything outside this run,
divide by the alignment length (count the bases per sequence in
results/panaroo/core_gene_alignment.aln) — the same 6,000 differences mean something very
different over a 0.5 Mb core than over a 2 Mb core.

**Relation to the panaroo section above.** The alignment length behind these counts is set by how many genes are core in this
particular sample set. The Panaroo paper shows how much that can move: on the same 319
Klebsiella pneumoniae genomes, Panaroo called 3372 core genes where Roary called 1800, and
for a clonal Mycobacterium tuberculosis outbreak the expected core was "approximately 4000"
genes. Read the SNP counts together with the Core/Soft core/Shell/Cloud table.

**"First 50 rows." note.** With a large set the report truncates the table to the first rows of the file (columns are
never cut). The complete matrix is always written to results/snp-dists/snp-dists.tsv — use
that file, not the on-screen table, for anything quantitative.

### Caveats

- These are raw counts with no evolutionary model behind them — no correction for
  transition/transversion bias and no correction for multiple substitutions at the same site
  — so treat the number as a ballpark similarity measure, not an evolutionary distance (a
  methodological caveat: snp-dists has no publication and neither paper makes this claim).

- Recombination is not accounted for either: a single imported tract can contribute hundreds
  of differences at once, so a large count can mean one horizontal transfer event rather
  than long independent divergence (again a methodological caveat, not a claim from these
  papers).

- The count depends entirely on how large the core genome is, which depends on which genomes
  you put in the run — adding one distant, contaminated or fragmented assembly shrinks the
  shared core for everybody and lowers all the counts, so values are not comparable between
  CompareM2 runs with different inputs; the Panaroo paper shows the same genomes yielding
  3372 versus 1800 core genes depending only on the tool used.

- Only the core is measured, so everything in the accessory genome is invisible here: two
  genomes can show 0 SNPs and still differ by a plasmid, a resistance cassette or dozens of
  genes — check the panaroo presence matrix and the AMRFinder section before calling two
  isolates the same.

- CompareM2 runs Panaroo with --clean-mode strict, which the paper describes as taking "a
  more aggressive approach to contamination and erroneous annotation removal" and which the
  authors note "can occasionally lead to rare plasmids being removed"; this shapes the gene
  set the alignment is built from.

### Quantitative record (8 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `9` | In the highly clonal M. tuberculosis London outbreak that Panaroo used as a control (413 genomes), the largest pairwise SNP distance in the whole dataset was 9 — an illustration of the scale of SNP counts within a single recent transmission chain. | “the maximum pairwise SNP distance within this dataset was 9” |
| `5 SNPs` | In that same outbreak, all but one isolate sat within 5 SNPs of the major clone. | “only one isolate being more than 5 SNPs from this” |
| `3372 and 3376` | The size of the core genome — the denominator behind every SNP count — depends heavily on the pangenome tool and its settings: Panaroo called 3372 (default) and 3376 (sensitive) core genes on the K. pneumoniae collection. | “sensitive modes, 3372 and 3376 respectively.” |
| `1800` | On exactly the same K. pneumoniae genomes, Roary's core genome was smaller by nearly a factor of two, showing how much the core alignment can shrink or grow without the genomes changing at all. | “Roary identified the smallest core genome of 1800 genes.” |
| `99%` | "Core" in that comparison meant a gene present in 99% of genomes — i.e. core is a presence threshold, not an absolute property of a gene. | “using the 99% presence threshold for core genes as used in Roary” |
| `approximately 4000` | For a clonal M. tuberculosis dataset the expected core genome was on the order of four thousand genes — a rough sense of how many genes typically sit behind a core-alignment SNP count. | “we would expect a very limited accessory genome and a core genome of approximately 4000” |
| `98%` | Genes are first grouped by CD-HIT at 98% nucleotide identity before any of Panaroo's graph corrections, so the sequences snp-dists compares within a column come from clusters built at that identity level. | “Panaroo runs CD-HIT (v4.8.1) at a high sequence identity threshold (98%)” |
| `smaller core genome estimates from unannotated genes` | Assembly fragmentation and missed annotations shrink the estimated core genome, and therefore the alignment these SNPs are counted over — one poor assembly lowers the counts for the entire set. | “genes being left unannotated resulting in smaller estimates of the core genome” |

### Not established from the papers (5)

- snp-dists has no peer-reviewed paper, so there is no published description of its
  algorithm, its defaults, or any validation of it — everything about the tool itself here
  comes from the command line CompareM2 runs, not from literature.

- Neither paper gives a SNP-distance threshold for 'same strain', 'same outbreak' or 'same
  clone'. The Panaroo paper reports one clonal M. tuberculosis outbreak with a maximum
  pairwise distance of 9, but it does not say those distances were computed from a core gene
  alignment, and outbreak cutoffs are species-specific and set by the epidemiological
  literature, not by this pipeline.

- How snp-dists treats gaps, Ns and other ambiguity codes at default settings is not
  documented in either paper; CompareM2 passes no flags, so whatever the installed version's
  default is applies, and I did not verify it against a running binary.

- The presence fraction a gene must reach to enter core_gene_alignment.aln is Panaroo's own
  default — CompareM2 does not set --core_threshold — and the paper states only the Roary-
  style 99% threshold it used for its own figures. So the exact core definition behind the
  alignment measured here is not established from the paper.

- The 'Core (≥99%)' row in the panaroo section of this report is computed by the report from
  the presence/absence matrix; it need not be the same gene set that Panaroo put into the
  alignment snp-dists read.

---

## fasttree
*Price MN, Dehal PS, Arkin AP 2010, PLoS ONE 5(3):e9490, doi:10.1371/journal.pone.0009490*

FastTree builds a phylogenetic tree from the core-gene alignment that Panaroo produced, so
you can see how your genomes are related by descent rather than just how similar they look
overall. It is designed to stay fast on large sets, and it buys that speed by searching tree
space less thoroughly than a full maximum-likelihood program.

**How it works.** FastTree 2 starts from a heuristic neighbour-joining tree, improves it with minimum-
evolution subtree-pruning-regrafting moves, and then rearranges branches under maximum
likelihood using nearest-neighbour interchanges only — it never does maximum-likelihood SPR
moves, which is why the authors call it "approximately-maximum-likelihood". Rate variation
across alignment columns is handled by the CAT approximation, which picks the single most
likely rate for each site out of 20 fixed values instead of integrating over a gamma
distribution; CompareM2 runs `FastTree -nt -gtr`, i.e. nucleotides under the general time-
reversible model, and FastTree's default SH-like local support values are computed at every
internal node.

### Reading the output

**The tree figure in the report.** CompareM2 draws the newick as a rectangular phylogram: each row is one genome, and
horizontal distance from left to right is accumulated branch length. Two things the figure
does not show — there is no scale bar, and internal node labels are not drawn. If you want
exact branch lengths or the support values, open fasttree/fasttree.newick.

**Branch length (horizontal distance).** Expected substitutions per site under GTR+CAT. Treat these as approximate: against branch
lengths re-optimised by PhyML under a 4-category gamma model, FastTree's internal branch
lengths correlated at r = 0.90, and for lengths between 0.01 and 1.0 the average percent
difference was 13%. Good enough to see which lineages are long-branched; not a calibrated
divergence estimate.

**Support value (internal node labels in fasttree.newick).** A number between 0 and 1 written immediately before the colon at each internal node, e.g.
`)0.987:0.00312`. It is an SH-like local support: FastTree takes the per-site likelihoods of
the current split and of the two nearest-neighbour-interchange alternatives around that
node, and applies the Shimodaira-Hasegawa test with 1,000 bootstrap replicates. It is a
local test of one node against two rivals, not a bootstrap of the whole alignment, so it
does not mean the same thing as a 0.99 bootstrap value.

**The leftmost node of the drawing.** An artefact of how newick is written, not a root. FastTree stores the tree with a
trifurcation at the root, and the paper states that the placement of the root is not
biologically meaningful and does not affect the likelihood. Do not read the leftmost split
as the deepest divergence; if you need direction, include an outgroup genome.

**Genome name colour and " · <cluster>" suffix.** These come from TreeCluster, if that step ran, not from FastTree. They colour the leaves by
cluster assignment and say nothing about how well supported the tree is.

### Caveats

- FastTree deliberately skips maximum-likelihood SPR moves, and the paper's own benchmark
  puts it below RAxML — 84.3% versus 88.4% of true splits on 5,000-sequence protein
  simulations — so a branching order you intend to build an argument on should be re-checked
  with a program that does a full SPR search.

- The support values are not bootstrap values: each one compares a node only against its two
  nearest-neighbour alternatives, and the authors themselves write that "local support
  values may be biased upwards because they do not consider all of" the alternate
  topologies, so a 1.00 here is weaker evidence than a 100% bootstrap.

- The paper explicitly warns that where nearby nodes are poorly resolved the supports should
  be interpreted cautiously, because a high-likelihood alternate topology may never have
  been considered — so a well-supported node sitting inside a bush of short branches is not
  trustworthy.

- Branch lengths come from the CAT approximation, and the paper states that if accurate
  branch lengths are essential then neither CAT nor the standard 4-category gamma model is
  sufficient; use them for topology and relative divergence, not as a molecular clock.

- The input is Panaroo's concatenated core-gene alignment, so the tree describes only the
  genes shared across the set and assumes one shared history for all of them — recombination
  and horizontal transfer break that assumption, and the paper does not address either
  (general methodological caution, not a FastTree finding).

### Quantitative record (13 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `84.3% → 85.0%` | Proportion of true splits recovered on simulated protein alignments with 5,000 sequences, for FastTree 2 as normally run and with its search heuristics removed | “alignments with 5,000 protein sequences from 84.3% to 85.0%.” |
| `88.4%` | RAxML 7.2.1 (which does maximum-likelihood SPR moves) on the same 5,000-sequence simulations | “was 88.4% accurate).” |
| `16%` | Of the true splits that PhyML 3 with SPR moves found but FastTree 2 missed, the fraction that were strongly supported (defined as SH-like and aLRT supports both at or above 95%) | “ratio test (aLRT) supports [21] of 95% or higher. Only 16% of the” |
| `20%` | Of the strongly supported splits found by PhyML 3 with SPRs but not by FastTree, the fraction that were actually incorrect | “but not FastTree, 20% were incorrect. Thus, few of the additional” |
| `96–98%` | Fraction of RAxML's splits with global bootstrap of 90% or higher that FastTree also recovered, on large genuine 16S alignments | “example, FastTree found 96–98% of RAxML” |
| `1,000` | Number of resamples used by the Shimodaira-Hasegawa test behind each local support value | “with 1,000 bootstrap replicates to estimate the confidence in the” |
| `r = 0.90` | Correlation between FastTree's SH-like local supports and PhyML's, on 250-sequence protein simulations | “between FastTree and PhyML (r = 0.90). For splits with local” |
| `0.008` | Average absolute difference between FastTree and PhyML support values, for splits supported at 0.9 or above by either tool | “average absolute difference was just 0.008, which is not much” |
| `0.880 instead of 0.887` | Ability of the support values to separate correct from incorrect splits, area under the ROC curve, FastTree versus PhyML | “the receiver operating curve (AOC) was 0.880 instead of 0.887” |
| `13%` | Average percent difference between FastTree's CAT-based internal branch lengths and PhyML's gamma-4 re-optimised lengths, for lengths between 0.01 and 1.0 | “average percent difference was 13%. For internal branch lengths” |
| `0.05 to 20` | Range of the 20 fixed per-site relative rates used by the CAT approximation | “the relative rates range from 0.05 to 20.” |
| `at least 100 times faster` | Speed of FastTree 2 relative to PhyML 3.0 and RAxML 7.2.1 on alignments of 500 sequences or more | “FastTree 2 is at least 100 times faster than either PhyML 3.0 or” |
| `76.2% → 78.0%` | Effect of the CAT rate model on accuracy for very small inputs — simulated protein alignments with just 10 sequences | “adding the CAT model improves FastTree’s accuracy from 76.2%” |

### Not established from the papers (5)

- The paper gives no recommended cutoff for calling a FastTree support value 'good'. The
  thresholds that appear in it (SH-like plus aLRT both at 95% or higher to label a PhyML
  split 'strongly supported'; bootstrap at 90% or higher for RAxML) are definitions used in
  its own comparisons, not advice for reading a FastTree tree.

- Every benchmark in the paper uses protein families (COGs) or 16S rRNA alignments. There is
  no benchmark on a concatenated bacterial core-gene alignment of the kind CompareM2 feeds
  it, and no statement about a minimum number of genomes or alignment length below which the
  approximations stop being adequate.

- The paper benchmarks FastTree 2.0.0 ("We used FastTree 2.0.0."); bioconda ships the 2.1.x
  line, whose changes the paper only sketches (better scaling, a -gamma option). The
  accuracy and speed figures above are therefore not measured on the exact binary CompareM2
  installs.

- catalogue.py declares threads=4 for this tool, but the command line invokes the plain
  `FastTree` binary rather than the OpenMP `FastTreeMP` build. The paper says only that
  FastTree 2.1 added parallel execution of the neighbour-joining phase, so whether those
  four threads are used at all is unverified.

- CompareM2 does not pass -gamma, -boot or any support-related flag, so the newick carries
  FastTree's default SH-like local supports. I did not verify against a real run that the
  values are in fact written for every internal node of a small (3-5 genome) tree.

---

## carveme
*Machado D, Andrejev S, Tramontano M, Patil KR 2018, Nucleic Acids Research 46:7542–7553, doi:10.1093/nar/gky537 (DIAMOND: Buchfink B, Reuter K, Drost H-G 2021, Nature Methods 18:366–368, doi:10.1038/s41592-021-01101-x)*

CarveMe turns each genome's predicted protein set into a genome-scale metabolic model: a
machine-readable list of the biochemical reactions that organism can probably run, wired
together into a network that a flux simulator can solve. It answers "what metabolic
capabilities does the annotation of this genome imply?" rather than "what does this organism
actually do in the lab".

**How it works.** CarveMe starts from one manually curated "universal" bacterial model assembled from the BiGG
database (elementally balanced, no blocked reactions, with a biomass equation), aligns your
proteins against the 30,814 BiGG-derived protein sequences with DIAMOND, converts alignment
scores into per-reaction confidence scores through gene-protein-reaction rules, and then
"carves": a mixed-integer linear program keeps high-scoring reactions, drops low-scoring
ones, and enforces network connectivity so no gappy dead ends remain. CompareM2 runs `carve
<bakta .faa> --output <sample>.xml` with no growth medium and no gap-filling, so only this
carving step happens; the optional medium-driven gap-filling, ensemble generation,
specialised Gram/archaea templates, experimental constraints and community merging are all
left at their defaults or unused.

### Reading the output

**Reactions.** Counted from every <reaction> element in the SBML, so it includes transporters, exchange
reactions, spontaneous reactions and the biomass equation — not just gene-associated ones
(the paper's own per-organism figures count "only gene-associated reactions", so your number
will be larger than the paper's). For scale, across CarveMe's 5587 RefSeq bacterial models
the average organism had 1308 reactions, the smallest was 238 (Mycoplasma ovis str Michigan)
and the largest 2472 (Klebsiella oxytoca str CAV1374). A few hundred reactions is expected
for a genuinely reduced genome; a few hundred for a free-living organism points at a poor or
truncated annotation.

**Metabolites.** Counted from <species> elements, which in SBML are compartment-specific: the same compound
in cytosol, periplasm and extracellular space counts three times. The curated universal
model has "2383 metabolites (representing 1503 unique compounds)", i.e. roughly one and a
half species entries per distinct chemical, so do not compare this column against the
paper's per-organism metabolite figures (average 792), which are counts of unique compounds.

**Genes.** The <fbc:geneProduct> entries — the genes from your own annotation that DIAMOND matched to a
BiGG gene and that survived carving. Across the 5587-model collection the average organism
had 691 metabolic genes, so this is a small subset of a bacterial CDS complement; the rest
of the genome simply has no counterpart in BiGG and is invisible to the model.

**Differences between genomes.** Treat between-genome differences in these three columns as annotation-plus-database
differences until shown otherwise. The paper states outright of its own model-size
statistics: "that these numbers are indicative, as they can be biased or restricted by the
quality of the gene annotation and the scope of the reaction database". Two strains
annotated with different tools, or one with a fragmented assembly, can differ by dozens of
reactions without any biological difference.

**"Not gap-filled for a specific medium".** The model is structurally simulation-ready — the carving problem carries a minimum growth-
rate constraint (default 0.1 h-1), so a carved model is built to be capable of producing
biomass — but nothing here was tuned to make it grow on a medium you care about. In the
paper's benchmark, four out of five test organisms did reproduce growth on minimal medium
from genome data alone with no gap-filling, and the fifth (Ralstonia solanacearum) needed
one added reaction (asparagine synthetase). So: use these models for network-level
comparison, presence/absence of pathways and as starting points for FBA; do not report a
predicted growth rate, auxotrophy or cross-feeding interaction from them without first
checking the model grows on a defined medium and curating it if it does not.

### Caveats

- The authors call these drafts themselves: the universal model is pre-curated and "can be
  readily used for simulation. Nonetheless, these models should still be considered as
  drafts subject to further refinement, as they might require organism-specific curation to
  reproduce certain phenotypes".

- Coverage is bounded by BiGG: the paper reports that primary metabolism is "essentially
  complete" while "peripheral pathways associated with secondary metabolism contain multiple
  gaps", so absence of a secondary-metabolite pathway in the model is not evidence the
  organism lacks it.

- Transport reactions are the known weak point — the paper attributes CarveMe's poorer
  substrate-utilisation performance for Bacillus subtilis to "the lack of annotated
  transporters", so uptake and secretion predictions (and hence cross-feeding claims) are
  the least trustworthy part of these models.

- CompareM2 passes no template flag, so the generic core biomass composition is used even
  for archaea and Gram-positives; the paper shows gene-essentiality sensitivity improves
  when Gram-specific or organism-specific biomass compositions are used instead, and
  separately provides an archaeal template "that contains ether lipids in the cell membrane
  composition and lacks peptidoglycans" which is not selected here.

### Quantitative record (13 verified)

| Value | Claim | Verbatim quote |
|---|---|---|
| `2383 metabolites (1503 unique compounds), 4383 reactions` | Size of the curated universal bacterial model that every carved model is cut out of | “2383 metabolites (representing 1503 unique compounds) and 4383 reactions” |
| `30 814` | Size of the reference protein database your proteins are aligned against with DIAMOND | “generate a sequence database with a total of 30 814 unique” |
| `691 metabolic genes, 1308 reactions, 792 metabolites` | Typical model size across the 5587 RefSeq bacterial models CarveMe reconstructed — the reference range for judging your own counts | “the average organism containing 691 metabolic genes, 1308 reactions and 792 metabolites” |
| `238 reactions` | Smallest network in that collection (a host-restricted genus), i.e. the low end of plausible reaction counts | “as few as 238 reactions (Mycoplasma” |
| `2472 reactions` | Largest network in that collection, i.e. the high end of plausible reaction counts | “2472 reactions (Klebsiella oxytoca str CAV1374)” |
| `four out of five` | Number of benchmark organisms whose ungapfilled model reproduced growth on minimal medium from genome data alone | “Notably, four out of five models generated with CarveMe” |
| `one single gap-filling reaction` | What the one failing benchmark model needed before it grew on minimal medium | “(R. solanacearum) required one single gap-filling reaction” |
| `0.1` | Default minimum growth rate (in h-1) enforced as a constraint inside the carving MILP | “growth rate (default: 0.1 h” |
| `30%` | How rare most reactions are across bacteria — relevant when a reaction is present in one of your genomes and absent in another | “with 30% of all reactions being present in less than 10%” |
| `r = –0.29, P = 0.0055` | Gap-filling need correlates negatively with genome size, which the authors attribute to poor annotation of transport genes in small genomes | “r = –0.29, P = 0.0055” |
| `3 min` | Wall-clock cost of one reconstruction as reported by the authors | “A single reconstruction required, on average, 3 min on a” |
| `~8,000-fold` | DIAMOND speedup over BLASTP in its least sensitive mode — the speed that makes per-genome carving cheap | “BLASTP (v2.10.0) we observed an ~8,000-fold speedup when” |
| `80-fold` | DIAMOND speedup when run at a sensitivity level matching BLASTP (ultra-sensitive), i.e. the cost of maximum sensitivity | “using the least sensitive mode, and still an 80-fold speedup when” |

### Not established from the papers (5)

- Which DIAMOND sensitivity mode carve invokes, and what identity/coverage/e-value cutoffs
  it applies to the alignments — neither paper states this, and the 2021 DIAMOND paper shows
  the four modes differ substantially in recall, so I cannot tell the reader how permissive
  the homology step is.

- What environment the carving MILP assumes when no medium is given. The paper states the
  MILP enforces a minimum growth rate (default 0.1 h-1) but does not say which uptake
  reactions are open during carving, so I cannot promise from the papers that a CompareM2
  model will carry biomass flux under any particular simulated medium — that has to be
  tested on the file.

- All MILP problems in the paper were solved with CPLEX 12.7, whereas the CompareM2
  environment solves with an open-source solver (per the note in catalogue.py). The paper
  says nothing about whether the carved network is solver-dependent when alternative optima
  exist; the paper does show that alternative equally plausible networks exist (that is why
  it offers ensembles), so this is worth checking rather than assuming.

- How many of a genome's CDS typically map into the model. The 691-gene average is the count
  of metabolic genes retained, not a fraction of the annotation, so the papers do not
  license a statement like "about X% of your genes are covered".

- Whether annotating with bakta versus prokka changes the resulting model. Not addressed by
  either paper.
