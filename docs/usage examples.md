
## Usage examples

Make a directory with the prokaryotic genomic assembly-files -or metagenomic bins- you want to investigate with Assemblycomparator2. 
Go into that directory in the terminal, and run the command `asscom2`. 
Assemblycomparator2 will then create a sub-directory, named "results_ac2/" containing a plethora of analysis results. 
  
  - Execute a "dry run". That is, to show what will be run without actually doing it.

    ```bash
    asscom2 --dry-run
    ```

  - Run Assemblycomparator2 on the genomes in the current directory:

    ```bash
    asscom2
    ```
    

##### A bit more advanced controls 

  - Run analyses that are relevant to metagenomic assemblies only (as opposed to isolates):

    ```bash
    asscom2 --until meta
    ```
    
  - Execute all jobs until one or more specific rules: (until implies including)
    
    ```bash
    asscom2 --until roary abricate
    ```
    
  - Select a specific MLST-scheme to use on all of the samples: (defaults to automatic)
    
    ```bash
    asscom2 --config mlst_scheme=hpylori
    ```
    
  - Select a specific roary blastp-identity: (default is 95)

    ```bash
    asscom2 --config roary_blastp_identity=90
    ```
      

