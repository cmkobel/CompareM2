# Usage


Overall, CompareM2 follows standard command line practices.
CompareM2 is built on top of Snakemake. Hence, when tweaking your run, you must pass the parameters through the `--config` key. Although all [Snakemake options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) are available to use, here we bring the ones that are necessary and useful for daily-driving CompareM2.


```txt
comparem2 [ --config KEY=VALUE [KEY2=VALUE]... ]
  [ --until RULE [RULE2]... ]
  [ --forcerun RULE [RULE2]... ]
  [ --downloads ]
  [ --printshellcmds ]
  [ --dry-run ]
  [ --status ]
  [ --version ]  [ --help ]  [ --cite]
```

## Usage examples

Here we bring som usage examples showcasing some of the more often used features.

  - Run *all* analyses across all fasta genome files in the current working directory.
    
    ```
    comparem2
    ```
    

  - Run only jobs *until* bakta
    
    ```
    comparem2 --until bakta
    ```

  - Run *all* analyses with specified input and output.
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" output_directory="my_analysis"
    ```

  - Use a *fofn* - a file of file names. 
    
    ```
    ls path/to/*.fna > my_fofn.txt; comparem2 --config fofn="my_fofn.txt"
    ```

  - Run a *dry run*.
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" --dry-run
    ```
  
  - Analyze some RefSeq references.
  
    ```
    comparem2 --config add_refseq="GCF_009734005.1,GCF_029023785.1"
    ```

  - Specify annotator. (default is "prokka")
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta"
    ```

  - Run only the *fast* rules. [(read more about pseudo rules)](https://comparem2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#pseudo-rules)
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta" --until fast
    ```

  - Run panaroo as well.
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta" --until fast panaroo
    ```

  - And pass a command line argument directly to panaroo.
    
    ```
    comparem2 --config input_genomes="path/to/genomes_*.fna" set_panaroo--threshold=0.95 annotator="bakta" --until fast panaroo 
    ```



## Options 

###  `--config KEY="VALUE" [KEY2="VALUE"]...`

The CompareM2 pipeline uses a configuration file with default settings to control the way it runs. This configuration can be changed on the fly -- This is to give the user easy access to the many features inside CompareM2. Multiple configuration options can be set at the same time as long as the all are stated after the `--config` statement on the command-line. The full default configuration file is available in ./config/config.yaml (in the installation directory). Below, we will document and showcase examples of the configuration parameters that are relevant when using CompareM2.

#### Configuration of input genomes


Input genomes for CompareM2 can be specified in several ways. The default is to use the `input_genomes` config key, and its default value is "*.fna *.fa *.fasta *.fas". In essence, this means that if the user does not specify this key, all genome fastas in the current working directory will be input into CompareM2 for analysis. The string given to input genomes is evaluated using GNU [LS(1)](https://www.linux.org/docs/man1/ls.html) which means that asterisks can be used for globbing. Examples below:
    
  - `input_genomes="*.fna *.fa *.fasta *.fas"` (default)
  - `input_genomes="path/to/my/genomes*.fna"`
  - `input_genomes="path/genome1.fna path/genome2.fna"`
  

##### File of file names
  
When analyzing larger sets of microbial genomes it can be useful to define these in a "file of file names" (fofn). This is supported by the `fofn` config key in CompareM2. When the `fofn` key is set, it always overrides the `input_genomes` key. A fofn-file can be generated simply by piping a list of filenames on the command line. Example below:


```bash

ls *.fna > fofn.txt
comparem2 --config fofn="fofn.txt"

```

  
##### Pre-annotated reference genomes

For many use cases it may be useful to add reference genomes from the [RefSeq or Genbank databases](https://www.ncbi.nlm.nih.gov/datasets/genome/), which contains consistently pre-annotated high quality genomes. Using the `add_refseq` config key, you can add one or more (comma separated) RefSeq or GenBank genome/assembly accessions as input to your CompareM2 run. The genomes and their [PGAP](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/process/) annotations will be automatically downloaded using the [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) command-line tools. The annotation is used as-is in the downstream tools in the CompareM2 pipeline. Examples below:

  - `add_refseq=<Refseq-Accession>[,<Refseq-Accession2>]...` 
  - `add_refseq="GCF_029023785.1"`
  
#### Configuration of the output directory

All results from CompareM2 are output into a subdirectory. The default name of this directory is "results_comparem2", but it can be changed, using the `output_directory` config key. Example below: 

  - `output_directory="results_comparem2"`.

#### Configuration of annotation tool

CompareM2 ships with two annotation tools: Bakta and Prokka. Bakta is default, and Prokka is bundled to better support users who work with Archaea. The choice of annotation tool in CompareM2 can have a large effect on your downstream results, as many tools use the output of the selected annotator as input. If you wish to change the annotation tool, you can do so with the `annotator` config key. Example below: 

  - `annotator="bakta"` bakta (default) | prokka
  

---

#### Passthrough arguments

From v2.8.2, CompareM2 has the ability to pass any command line argument (option-parameter pair) through to any rule in the workflow. This is done by using a generalized "passthrough argument" feature that recognizes config argument options starting with string "set_" and passes them to the correct rule upon generating the shell scripts for each rule in the workflow. The general syntax for these passthrough arguments is `set_<rule><option>=<parameter>` where rule is the rule name, option is the option key, and parameter is the parameter value. 

!!! note
    This feature requires modification of Snakemake such that it can accept special characters through the config strings given at the command line. This modification can easily be done using the following command that ships with the bioconda package:
    ```
    enable_passthrough_parameters_comparem2
    ```
    
    Otherwise you might receive the Snakemake error: "Invalid config definition: Config entry must start with a valid identifier."

An example can be used to explain how this feature can be used in practice: Consider using the Prokka annotator, which is capable of annotating both bacterial and archaeal genomes. By default, Prokka is set to bacterial annotation, so in case we want to annotate an archaea, we can set the "--kingdom" argument to "archaea". In this case the rule name is `prokka`, the option key is `--kingdom` and the parameter value is `archaea`. When using CompareM2, this setting can be set following the passthrough argument syntax like so:
 
```bash
# comparem2 --until set_<rule><key>=<value> # Syntax template.
comparem2 --config set_prokka--kingdom=archaea
```

Notice how the double dash prefix in "--kingdom" is part of the the set_ string. This is because many different styles of command line argument options need to be supported (e.g.: "--command_key", "--command-key", "-command_key" etc). 

In some cases, command line options are flags, meaning that they need no parameter value. In this case, an empty string can be given as parameter value:

```bash
comparem2 --config set_prokka--rfam="" # --rfam enables searching for ncRNAs with Infernal+Rfam.
```

In case of non-empty parameter values, use of apostrophes is optional.

Using a space separator, several command line arguments can be given at once for several different tools. In the following example we're also loosening the Panaroo core genome identity "--threshold" option down to 95% to increase the apparent number of genes in the core genome.

```bash
comparem2 --config set_prokka--kingdom=archaea set_panaroo--threshold=0.95 --until panaroo fast
```

  
CompareM2 comes with a number of sane default arguments which can be observed [here](https://github.com/cmkobel/CompareM2/blob/master/config/config.yaml). Any passthrough argument that the user gives on the command line overwrites these defaults.


##### Validating command line arguments

There are no limitations on which command line arguments can be passed to the passthrough argument feature. Thus, when modifiying the options using the passthrough arguments feature, the user should follow the documentation of each individual tool to make sure that the command line arguments given are valid. In order to validate that the arguments given to rules are as expected, the full generated shell command of each rule can be printed with `-p`. It is especially useful to do this in conjunction with the `--dry-run` argument. Example with Panaroo below:

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




---

### General Snakemake commands

As CompareM2 is built on top of Snakemake, it supports its built in features and command line arguments. Below is a listing of the most relevant Snakemake commands to use together with CompareM2.
            
#### `--until RULE [RULE2]...`
Select to run up until and including a specific rule in the rule graph. Available rules:
abricate annotate antismash assembly_stats bakta checkm2 copy dbcan eggnog fasttree gapseq gapseq_find gtdbtk interproscan iqtree kegg_pathway mashtree mlst prokka sequence_lengths snp_dists treecluster antismash_download bakta_download checkm2_download dbcan_download eggnog_download gtdb_download panaroo
    
There are also a number of pseudo rules, effectively "shortcuts" to a list of rules.
  - downloads   (Run rules that download and setup up necessary databases.)
  - fast        (Only rules that complete within a few seconds. Useful for testing.)
  - isolate     (Only rules that are relevant for genomes of isolate origin.)
  - meta        (Only rules that are relevant for genomes "MAGs" of metagenomic origin.)
  - report      (Re-renders the report.)

---
          
#### `--forcerun RULE [RULE2]...`

Force rerunning of one or more rules that already have been completed. This is generally necessary when changing running parameters in the config (see "--config" above).

---

#### `--downloads`

Download all databases without performing any analyses.

---

#### `--printshellcmds`, `-p`

Print the full generated shell commands of each rule in the workflow. 

---
    
#### `--dry-run`

Run a "dry run": Shows what will run without doing it.

---

#### `--status`

Print the state of completion of the rules in the pipeline. The percentage of completed files are shown.

Run `comparem2 --status` at any time to get an overview of what has been done and what is missing in any project directory. Should be run from the same directory as where comparem2 was run. 

---

#### `--version`, `-v `

Show current version.

---
    
#### `--help`, `-h`

Show this help and exit.

---
 
#### `--cite`

Show citation info.

---
        
## Environment variables
No environment variables are strictly necessary to set, but the following might be useful:

  - `COMPAREM2_PROFILE` (default "profile/apptainer/local") specifies which of the Snakemake profiles to use. This can be useful for running CompareM2 on a HPC or using specific settings on a large workstation. Check out the bundled profiles in path profile/* (possibly in $CONDA_PREFIX/comparem2/profile/\*).
  
  - `COMPAREM2_DATABASES` (default "databases/") specifies a database location. Useful when sharing a database installation between various users on the same workstation or HPC.



## Output
CompareM2 creates a directory named "results_comparem2/" (or what the output_directory parameter is set to) that contains all of the analysis results that are computed.

Results from input genomes are in dir "samples/" and results across all samples are in the root. 

The report is named "report_&lt;title>.html" after the title of the run which defaults to the name of the current working directory.



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
│     ├── antismash/
│     ├── bakta/
│     ├── dbcan/
│     ├── eggnog/
│     ├── <sample>.fna
│     ├── interproscan/
│     ├── prokka/
│     └── sequence_lengths/
├── snp-dists/
├── tables/
├── treecluster/
└── version_info.txt
```

For the file tree of each of the analysis tools, please consult the respective documentation.

{!resources/footer.md!}
