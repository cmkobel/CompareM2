# Usage

## Usage examples

Overall, Assemblycomparator follows standard command line practices.
Assemblycomparator2 is built on top of Snakemake. Hence, when customizing the Assemblycomparator pipeline you must pass the parameters through the `--config` key. All [Snakemake options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) are freely available for advanced users.


  - Run *all* analyses with specified input genomes.
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna"
    ```

  - Use a *fofn* - a file of file names. 
    
    ```
    asscom2 --config fofn="my_fofn.txt"
    ```

  - Use custom output dir. (default is "results_ac2")
    
    ```
    asscom2 --config input_genomes="path/to/genomes_*.fna" output_directory="my_analysis"
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




## Full `--help`

```txt
{!docs/help_text.txt!}
```


## Demo reports

These demo reports are available for inspiration.

  - [report_strachan_campylo.html](https://github.com/cmkobel/assemblycomparator2/raw/master/tests/strachan_campylo/report_strachan_campylo.html.zip)

    32 Campylobacter genomes, Metagenome and genome sequencing from the rumen epithelial wall of dairy cattle. From Nature 2022 - Strachan et al. (doi.<nolink />org/10.1038/s41564-022-01300-y).
  - 





