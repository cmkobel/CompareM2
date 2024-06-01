# Usage

## Usage examples

Make a directory with the prokaryotic genomic assembly-files -or metagenomic bins- you want to investigate with Assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
Assemblycomparator2 will then create a sub-directory, named "results_ac2/" containing a plethora of analysis results. 
  
  - Execute a "dry run". That is, to show what will be run without actually doing it.

    ```
    asscom2 --dry-run
    ```

  - Run Assemblycomparator2 on the genomes in the current directory:

    ```
    asscom2
    ```
    

## A bit more advanced controls 

  - Run analyses that are relevant to metagenomic assemblies only (as opposed to isolates):

    ```
    asscom2 --until meta
    ```
    
  - Execute all jobs until one or more specific rules: (until implies including)
    
    ```
    asscom2 --until roary abricate
    ```
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    ```
    asscom2 --config mlst_scheme=hpylori
    ```
    
  - Select a specific roary blastp-identity: (default is 95)

    ```
    asscom2 --config roary_blastp_identity=90
    ```
      



## Full `--help`

```txt
{!docs/help_text.txt!}
```


## Demo reports

These demo reports are available for inspiration.

  - [report_strachan_campylo.html](https://github.com/cmkobel/assemblycomparator2/raw/master/tests/strachan_campylo/report_strachan_campylo.html.zip)

    32 Campylobacter genomes, Metagenome and genome sequencing from the rumen epithelial wall of dairy cattle. From Nature 2022 - Strachan et al. (doi.<nolink />org/10.1038/s41564-022-01300-y).
  - 





