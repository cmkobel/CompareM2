# Usage


Overall, CompareM2 follows standard command line practices.
CompareM2 is built on top of Snakemake. Hence, when tweaking your run, you must pass the parameters through the `--config` key. All [Snakemake options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) are available as well.


```txt
comparem2 [ --config KEY=VALUE [KEY2=VALUE]... ]
  [ --until RULE [RULE2]... ]
  [ --forcerun RULE [RULE2]... ]
  [ --dry-run ]
  [ --version ]  [ --help ]
```

## Usage examples


  - Run *all* analyses across all fasta files in the current working directory.
    
    ```
    comparem2
    ```

  - Run only jobs *until* prokka
    
    ```
    comparem2 --until prokka
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

###  `--config KEY=VALUE [KEY2=VALUE]...`
Pass a parameter to the snakemake pipeline, where the following keys are available, defaults are stated as standard.
    
  - `input_genomes="*.fna *.fa *.fasta *.fas"` Path to input genomes. As the default value indicates, all fasta type files in the present directory will be analyzed.

  - `fofn="fofn.txt"` Deactivated by default. When set to a path it overrides key input_genomes. A fofn can be created with `ls *.fna > fofn.txt`

  - `output_directory="results_comparem2"` All results are written here.

  - `annotator="prokka"` Choice of annotation tool. Alternatively "bakta".
    



#### Passthrough arguments

From v2.8.2, CompareM2 has the ability to pass _any_ command line argument (option parameter pair) through to _any_ rule in the workflow. This is done by using a generalized "passthrough argument" feature that recognizes config arguments options starting with string "set_" and passes them to the correct rule upon generating the shell scripts for each rule in the workflow. The general syntax for these passthrough arguments is `set_<rule><option>=<parameter>` where rule is the rule name, option is the option key and parameter is the parameter value. 

!!! info
    This feature requires modification of Snakemake such that it can accept special characters through the config strings given at the command line. This modification can easily be done using the following command that ships with the bioconda package:
    ```
    enable_passthrough_parameters_comparem2
    ```

An example can be used to explain how this feature can be used in practice: Consider using the Prokka annotator (rule name `prokka`), which is capable of annotating both bacterial and archaeal genomes it is necessary to set the "--kingdom" argument to "archaea" when analyzing archaea. This is for Prokka to use the correct internal databases when annotating the genes, because the default is bacteria. In this case the rule name is `prokka`, the option key is `--kingdom` and the parameter value is `archaea`. When using CompareM2, this setting can be set following the passthrough parameter syntax like so:
 
```
# comparem2 --until set_<rule><key>=<value> # Syntax template.
comparem2 --until set_prokka--kingdom=archaea
```

Notice how the double dash prefix is part of the the set_ string. This is because many different styles of command line argument options need to be supported (e.g.: `--command_key`, --command-key, -command_key` etc).

In some cases, command line arguments are options that work like a flag, meaning that they need no parameter value. In this case, an empty string can be given an the parameter value:

```
comparem2 --until set_prokka--compliant=""
```

In case of non-empty parameter values, it is optional to use apostrophes.

Several command line arguments can be given at once for several different tools, using a space separator. In the following example we're also loosening the Panaroo core genome identity `--threshold` option down to 95% to increase the apparent number of genes in the core genome.

```
comparem2 --config set_prokka--kingdom=archaea set_panaroo--threshold=0.95 --until panaroo fast
```


!!! note 
    Remember that there are no limitations on which command line arguments can be passed to the passthrough argument feature. The user should follow the documentation of each individual tool to make sure that the command line arguments given are valid.

  
CompareM2 comes with a number of sane default arguments which can be observed [here](https://github.com/cmkobel/CompareM2/blob/master/config/config.yaml). Any passthrough argument that the user gives on the command line overwrites these defaults.


            
### `--until RULE [RULE2]...`
Select to run up until and including a specific rule in the rule graph. Available rules:
abricate annotate antismash assembly_stats bakta busco checkm2 copy dbcan eggnog fasttree gapseq gapseq_find gtdbtk interproscan iqtree kegg_pathway mashtree mlst prokka sequence_lengths snp_dists treecluster antismash_download bakta_download busco_download checkm2_download dbcan_download eggnog_download gtdb_download panaroo
    
There are also a number of pseudo rules, effectively "shortcuts" to a list of rules.
  - downloads   (Run rules that download and setup up necessary databases.)
  - fast        (Only rules that complete within a few seconds. Useful for testing.)
  - isolate     (Only rules that are relevant for genomes of isolate origin.)
  - meta        (Only rules that are relevant for genomes "MAGs" of metagenomic origin.)
  - report      (Re-renders the report.)
          
### `--forcerun RULE [RULE2]...`
Force rerunning of one or more rules that already have been completed. This is generally necessary when changing running parameters in the config (see "--config" above).
    
### `--dry-run`
Run a "dry run": Shows what will run without doing it.

### `--version`, `-v `
Show current version.
    
### `--help`, `-h`
Show this help and exit.
        
## Environment variables
No environment variables are strictly necessary to set, but the following might be useful:

  - `COMPAREM2_PROFILE` (default "profile/apptainer/local") specifies which of the Snakemake profiles to use. This can be useful for running CompareM2 on a HPC or using specific settings on a large workstation. Check out the bundled profiles in path profile/* (possibly in $CONDA_PREFIX/comparem2/profile/\*).
  
  - `COMPAREM2_DATABASES` (default "databases/") specifies a database location. Useful when sharing a database installation between various users on the same workstation or HPC.
  
## Output
Creates a directory named "results_comparem2/" (or what the output_directory parameter is set to) that contains all of the analysis results that are computed.

A file tree with depth level 1 looks like so:

```txt
results_comparem2/
├── abricate/
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
├── snp-dists/
├── tables/
├── treecluster/
└── version_info.txt
```

Results from input genomes are in dir "samples/" and results across all samples are in the root. The report is named after the title of the run which is the same as the name of the current working directory.

<sample\> is defined as the basename without extension of each file given as input genome.

The samples/<sample\> directory for each sample looks like so: (Again, depth level 1 only.)
```txt
results_comparem2/
└── samples/
   └── <sample>/
      ├── antismash/
      ├── bakta/
      ├── busco/
      ├── dbcan/
      ├── eggnog/
      ├── <sample>.fna
      ├── interproscan/
      ├── prokka/
      └── sequence_lengths/
```

For the file tree of each of the analysis tools you will have to consult the documentation of each tool.



{!resources/footer.md!}
