# Usage


Overall, Assemblycomparator follows standard command line practices.
Assemblycomparator2 is built on top of Snakemake. Hence, when tweaking your run, you must pass the parameters through the `--config` key. All [Snakemake options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) are available as well.


```txt
asscom2 [ --config KEY=VALUE [KEY2=VALUE]... ]
  [ --until RULE [RULE2]... ]
  [ --forcerun RULE [RULE2]... ]
  [ --dry-run ]
  [ --version ]  [ --help ]
```

## Usage examples




  - Run *all* analyses across all fasta files in the current working directory.
    
    ```
    asscom2
    ```

  - Run only jobs *until* prokka
    
    ```
    asscom2 --until prokka
    ```

  - Run *all* analyses with specified input and output.
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" output_directory="my_analysis"
    ```

  - Use a *fofn* - a file of file names. 
    
    ```
    asscom2 --config fofn="my_fofn.txt"
    ```

  - Run a *dry run*.
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" --dry-run
    ```

  - Specify annotator. (default is "prokka")
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta"
    ```

  - Run only the *fast* rules. [(read more about pseudo rules)](https://assemblycomparator2.readthedocs.io/en/latest/30%20what%20analyses%20does%20it%20do/#pseudo-rules)
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta" --until fast
    ```

  - Run panaroo as well.
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" annotator="bakta" --until fast panaroo
    ```


## Options 

###  `--config KEY=VALUE [KEY2=VALUE]...`
Pass a parameter to the snakemake pipeline, where the following keys are available, defaults are stated as standard. (Longer explanation in paranthesis.)
    
  - `input_genomes="*.fna *.fa *.fasta *.fas"` (Path to input genomes. As the default value indicates, all fasta type files in the present directory will be analyzed).

  - `fofn="fofn.txt"` (Deactivated by default. When set to a path it overrides key input_genomes.)

  - `output_directory="results_ac2"` (All results are written here.)

  - `annotator="prokka"` (Choice of annotation tool. Alternatively "bakta".)
    
There also is a series of parameters inside specific tools that can be tweaked directly from the command line interface. Please consult the individual documentation of the relevant tool for more info. Examples of possible values are given in parenthesis:

  - `mlst_scheme="automatic"`
  - `prokka_rfam=true` (true or false)
  - `prokka_compliant=true` (true or false)
  - `treecluster_threshold=0.045` (interpreted as float)
  - `iqtree_bootstraps=100` (interpreted as int)
            
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

  - `ASSCOM2_PROFILE` (default "profile/apptainer/local") specifies which of the Snakemake profiles to use. This can be useful for running Assemblycomparator2 on a HPC or using specific settings on a large workstation. Check out the bundled profiles in path profile/* (possibly in $CONDA_PREFIX/assemblycomparator2/profile/\*).
  
  - `ASSCOM2_DATABASES` (default "databases/") specifies a database location. Useful when sharing a database installation between various users on the same workstation or HPC.
  
## Output
Creates a directory named "results_ac2/" (or what the output_directory parameter is set to) that contains all of the analysis results that are computed.

A file tree with depth level 1 looks like so:

```txt
results_ac2/
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
├── report_strachan_campylo.html
├── samples/
├── snp-dists/
├── tables/
├── treecluster/
└── version_info.txt
```

Results from input genomes are in dir "sample/" and results across all samples are in the root. The report is named after the title of the run which is the same as the name of the current working directory.



{!resources/footer.md!}
