# Usage

CompareM2 is built on top of Snakemake. Pipeline parameters are passed via `--config`, and all [Snakemake command-line options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) are available.

```txt
comparem2 [ --config KEY=VALUE [KEY2=VALUE]... ]
  [ --until RULE [RULE2]... ]
  [ --forcerun RULE [RULE2]... ]
  [ --downloads ]
  [ --printshellcmds ]
  [ --dry-run ]
  [ --status ]
  [ --version ]  [ --help ]  [ --cite ]
```


## Examples

Run all analyses on all genome files (`*.fna *.fa *.fasta *.fas`) in the current directory:

```bash
comparem2
```

Run only the fast analyses:

```bash
comparem2 --until fast
```

Specify input and output paths:

```bash
comparem2 --config input_genomes="path/to/genomes_*.fna" output_directory="my_analysis"
```

Use a file-of-filenames (fofn):

```bash
ls path/to/*.fna > my_fofn.txt
comparem2 --config fofn="my_fofn.txt"
```

Dry run (preview what will run without executing):

```bash
comparem2 --config input_genomes="path/to/genomes_*.fna" --dry-run
```

Analyze NCBI reference genomes by accession:

```bash
comparem2 --config add_ncbi="GCF_009734005.1,GCF_029023785.1"
```

Use Prokka instead of the default Bakta annotator:

```bash
comparem2 --config annotator="prokka"
```

Combine options — run fast analyses plus panaroo with Bakta annotation:

```bash
comparem2 --config input_genomes="path/to/genomes_*.fna" --until fast panaroo
```

Pass a parameter directly to an underlying tool:

```bash
comparem2 --config set_panaroo--threshold=0.95 --until fast panaroo
```


## Configuration reference

### Input genomes (`input_genomes`)

A glob pattern specifying which genome files to analyze. The default picks up all common FASTA extensions in the current directory:

  - `input_genomes="*.fna *.fa *.fasta *.fas"` (default)
  - `input_genomes="path/to/my/genomes*.fna"`
  - `input_genomes="path/genome1.fna path/genome2.fna"`

### File of file names (`fofn`)

For larger sets of genomes, list paths in a text file (one per line). When set, `fofn` overrides `input_genomes`:

```bash
ls *.fna > fofn.txt
comparem2 --config fofn="fofn.txt"
```

### Pre-annotated NCBI reference genomes (`add_ncbi`)

Add reference genomes from [NCBI/GenBank](https://www.ncbi.nlm.nih.gov/datasets/genome/) by accession. Genomes and their [PGAP](https://www.ncbi.nlm.nih.gov/ncbi/annotation_prok/process/) annotations are downloaded automatically via [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/). Multiple accessions are comma-separated:

```bash
comparem2 --config add_ncbi="GCF_029023785.1,GCF_009734005.1"
```

### Output directory (`output_directory`)

Where results are written. Default: `results_comparem2`.

```bash
comparem2 --config output_directory="my_results"
```

### Annotation tool (`annotator`)

CompareM2 ships with two annotators:

  - **`bakta`** (default) — recommended for bacteria
  - **`prokka`** — also supports archaea

The choice of annotator affects many downstream analyses, as tools like panaroo, eggnog, and interproscan consume its output.

```bash
comparem2 --config annotator="prokka"
```

### Report title (`title`)

Custom title for the HTML report. Defaults to the name of the current working directory.


## Passthrough arguments

CompareM2 can forward arbitrary command-line arguments to any underlying tool using a `set_` prefix in the config. The syntax is:

```
set_<rule><option>=<value>
```

Where `<rule>` is the Snakemake rule name, `<option>` is the tool's command-line flag (including dashes), and `<value>` is the parameter value.

**Example:** Set Prokka's `--kingdom` flag to `archaea`:

```bash
comparem2 --config set_prokka--kingdom=archaea
```

**Flag-only arguments** (no value) use an empty string:

```bash
comparem2 --config set_prokka--rfam=""
```

**Multiple passthrough arguments** can be combined:

```bash
comparem2 --config set_prokka--kingdom=archaea set_panaroo--threshold=0.95 --until panaroo fast
```

The [default passthrough arguments](https://github.com/cmkobel/CompareM2/blob/master/config/config.yaml) can be overridden by specifying the same key on the command line.

!!! note
    Passthrough arguments require a modification to Snakemake to accept special characters in config strings. Run this command (included in the bioconda package) to enable it:
    ```
    enable_passthrough_parameters_comparem2
    ```
    Without this, Snakemake will report: "Invalid config definition: Config entry must start with a valid identifier."

### Validating passthrough arguments

Use `-p --dry-run` to preview the generated shell commands and verify that your arguments are being passed correctly:

```bash
comparem2 --config set_panaroo--threshold=0.99 --until panaroo -p --dry-run
#> [...]
#>    panaroo \
#>        -o results_comparem2/panaroo \
#>        -t 16 \
#>        --clean-mode sensitive \
#>        --core_threshold 0.95 \
#>        --threshold 0.99 \
#> [...]
```


## Snakemake options

### `--until RULE [RULE2]...`

Run only the specified rule(s) and their dependencies. Multiple rules can be listed. Available rules:

`abricate` `amrfinder` `annotate` `antismash` `assembly_stats` `bakta` `bootstrap_mashtree` `carveme` `checkm2` `copy` `dbcan` `eggnog` `fasttree` `gapseq_find` `gapseq_fill` `gtdbtk` `interproscan` `iqtree` `kegg_pathway` `mashtree` `mlst` `panaroo` `prokka` `sequence_lengths` `snp_dists` `treecluster`

Download rules: `antismash_download` `bakta_download` `checkm2_download` `dbcan_download` `eggnog_download` `gtdb_download`

### Pseudo-rules

Pseudo-rules are shortcuts that run a curated set of analyses:

| Pseudo-rule | Description | Included analyses |
|---|---|---|
| `fast` | Completes in seconds; useful for testing | sequence_lengths, assembly-stats, mashtree |
| `meta` | Analyses relevant for MAGs | annotation, assembly-stats, sequence_lengths, checkm2, eggnog, kegg_pathway, dbcan, interproscan, gtdbtk, mashtree |
| `isolate` | Analyses relevant for clinical isolates | annotation, assembly-stats, sequence_lengths, eggnog, kegg_pathway, gtdbtk, mlst, amrfinder, panaroo, fasttree, snp-dists, mashtree |
| `downloads` | Download and set up all databases | All database download rules |
| `report` | Re-render the report only | Report generation |

Usage: `comparem2 --until meta` or `comparem2 --until isolate`

### `--forcerun RULE [RULE2]...`

Force re-execution of completed rules. Necessary when changing config parameters for a rule that has already run.

### `--printshellcmds`, `-p`

Print the generated shell command for each rule.

### `--dry-run`

Show what would run without executing anything.


## CompareM2-specific options

These options do not invoke the Snakemake pipeline.

| Option | Description |
|---|---|
| `--downloads` | Download all databases without running analyses |
| `--status` | Show completion status of each rule in the current project directory |
| `--version`, `-v` | Show version |
| `--help`, `-h` | Show help |
| `--cite` | Show citation information |


## Output structure

CompareM2 writes all results to the output directory (default: `results_comparem2/`). Per-sample results are in `samples/<sample>/`, and cross-sample results are in the root.

The report is named `report_<title>.html`, where `<title>` defaults to the current working directory name.

```txt
results_comparem2/
├── amrfinder/
├── assembly-stats/
├── benchmarks/
├── checkm2/
├── fasttree/
├── gtdbtk/
├── iqtree/
├── kegg_pathway/
├── mashtree/
├── metadata.tsv
├── mlst/
├── panaroo/
├── report_<title>.html
├── samples/
│  └── <sample>/
│     ├── <sample>.fna
│     ├── antismash/
│     ├── bakta/
│     ├── dbcan/
│     ├── eggnog/
│     ├── gapseq/
│     ├── interproscan/
│     ├── prokka/
│     └── sequence_lengths/
├── snp-dists/
├── visuals/
├── treecluster/
└── version_info.txt
```

For per-tool file details, consult the respective tool's documentation.

{!resources/footer.md!}
