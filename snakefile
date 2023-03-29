

# snakemake --snakefile ~/assemblycomparator2/snakefile --profile ~/assemblycomparator2/configs/slurm/ --cluster-config ~/assemblycomparator2/configs/cluster.yaml 

__version__ = "v2.4.0" # Ttree places to bump. Here, in the bottom of the report, in the report snakefile
__author__ = 'Oliver Kjærlund Hansen & Carl M. Kobel'

import os
from os import listdir
from os.path import isfile, join
#import yaml
import pandas as pd
import numpy as np
from shutil import copyfile
#import time
#import re
#from shutil import copyfile
#import re

# For void_report
import subprocess
import datetime
        


# --- Read important variables -----------------------------------------------
cwd = os.getcwd()
batch_title = cwd.split("/")[-1]
base_variable = os.environ['ASSCOM2_BASE'] # rename to ASSCOM2_BASE

print("/*")
print()
print("         █████╗ ███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ ")
print("        ██╔══██╗██╔════╝██╔════╝██╔════╝██╔═══██╗████╗ ████║╚════██╗")
print("        ███████║███████╗███████╗██║     ██║   ██║██╔████╔██║ █████╔╝")
print("        ██╔══██║╚════██║╚════██║██║     ██║   ██║██║╚██╔╝██║██╔═══╝ ")
print("        ██║  ██║███████║███████║╚██████╗╚██████╔╝██║ ╚═╝ ██║███████╗")
print("        ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝")
print("                       A.K.A. assemblycomparator2                   ")
print("                         Please log issues at:                      ")
print("              github.com/cmkobel/assemblycomparator2/issues         ")
print("                                                                    ")
print(f"    batch_title:           {batch_title}")
print(f"    roary_blastp_identity: {config['roary_blastp_identity']} (default 95)")
print(f"    mlst_scheme:           {config['mlst_scheme']} (default automatic)   ")
print(f"    base_variable:         {base_variable}                               ")
#print(f"    kraken2 database:      {config['asscom2_kraken2_db']}                ")
#print(f"    gtdb:                  {config['gtdbtk_data_path']}                  ")


#results_directory = "output_asscom2"



results_directory = "results_ac2"



#reference = config["reference"]



# --- Read in relevant files in the current working directory ----------------

relative_wd = "."
#extension_whitelist = ["fna", "fa", "fas", "fasta", "seq"] # old 
extension_whitelist = ["fna", "fa", "fas", "fasta", "seq", "gb", "fq", "gff", "gfa", "clw", "sth", "gz", "bz2"]

present_files = [f for f in listdir(relative_wd) if isfile(join(relative_wd,f))]


df = pd.DataFrame(data = {'input_file': present_files})

# Check that the directory is not empty.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...")
    raise Exception("Zero genomic files present.")



df = df[~df["input_file"].str.startswith(".", na = False)] # Remove hidden files
df['sample_raw'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()] # Extract everything before the extension dot.
df['sample'] = df['sample_raw'].str.replace(' ','_').str.replace(',','_') # Convert spaces and commas to underscore 
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()] # Extract extension
df['input_file_fasta'] = results_directory + "/samples/" + df['sample'] + "/" + df['sample'] + ".fa" # This is where the input file is copied to in the first snakemake rule.

df = df[df['extension'].isin(extension_whitelist)] # Remove files with unsupported formats.

# Check that the directory is not empty, again.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...(2)")
    raise Exception("Zero genomic files present.")


#df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscores.

  
# --- Displaying filtered dataframe ready for analysis ------------------------
print() # Padding
df = df.reset_index(drop = True)
#print(df[['input_file', 'sample', 'extension']])
print(df[['input_file', 'extension', 'input_file_fasta']])
print("//")
print()




# --- Make sure the output directory exists. ----------------------------------

""" try: 
    os.mkdir(f"{results_directory}") # If running with local profile, the directory won't be created. I'm not sure if it needs to though?
except:
    pass
 """


# The modification time of this file tells the report subpipeline whether it needs to run. Thus, void_report is called in the end of every successful rule.
#void_report = f"touch {results_directory}/.asscom2_void_report.flag"
void_report = f"date -Iseconds >> {results_directory}/.asscom2_void_report.flag"






# --- Collect all targets. ------------------------------------------
rule all:
    input: expand([\
        "{results_directory}/metadata.tsv", \
        "{results_directory}/.install_report_environment_aot.flag", \
        "{results_directory}/checkm2/quality_report.tsv", \
        "{results_directory}/assembly-stats/assembly-stats.tsv", \
        "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", \
        "{results_directory}/roary/summary_statistics.txt", \
        "{results_directory}/abricate/card_detailed.tsv", \
        "{results_directory}/mashtree/mashtree.newick", \
        "{results_directory}/mlst/mlst.tsv", \
        "{results_directory}/fasttree/fasttree.newick", \
        "{results_directory}/gtdbtk/gtdbtk.bac.summary.tsv", \
        "{results_directory}/snp-dists/snp-dists.tsv", \
        "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"], \
        results_directory = results_directory, sample = df["sample"]) 



# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
rule copy:
    #input: "{sample}"
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: "{results_directory}/samples/{sample}/{sample}.fa"
    #log: "logs/{results_directory}_{wildcards.sample}.out.log"
    container: "docker://pvstodghill/any2fasta"
    conda: "conda_definitions/any2fasta.yaml"
    resources:
        mem_mb = 512,
        runtime = "00:10:00"
    shell: """

        any2fasta {input.genome:q} > {output}

    """  


# Write the df table to the directory for later reference.
# Why isn't this a run: instead of a shell: ?
rule metadata:
    input: expand("{results_directory}/samples/{sample}/{sample}.fa", results_directory = results_directory, sample = df["sample"]) # From rule copy
    output: "{results_directory}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    resources:
        runtime = "01:00:00"
    #run:
        #df.to_csv(str(output), index_label = "index", sep = "\t")
        #os.system(f"cp ${{ASSCOM2_BASE}}/scripts/{report_template_file_basename} {results_directory}")
    shell: """

        echo '''{params.dataframe}''' > {output}


        {void_report}
    """



# --- Downloads -----------------------------------------------------

# This rule runs once, downloading the busco dataset that is needed for rule busco_individual.
# Make sure that this job is run on a node that has internet access.
rule busco_download:
    output:
        #touch("{base_variable}/databases/busco/busco_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        touch("{base_variable}/databases/busco/file_versions.tsv") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
    conda: "conda_definitions/busco.yaml"
    shell: """

        
        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Checking for internet access using google."
            ping -q -c1 google.com &>/dev/null && echo "Online" || echo "Warning: It seems like you don't have internet access? Downloading will probably fail."

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            # Busco is a bit stupid in the way that it requires an input file, but doesn't read it when you just download.
            touch dummy.fasta
            
            # https://busco.ezlab.org/busco_userguide.html#download-and-automated-update
            busco \
                --in dummy.fasta \
                --out shouldnotbenecessarytosettheoutdirwhenjustdownloading \
                --mode geno \
                --auto-lineage-prok \
                --force \
                --download_path {wildcards.base_variable}/databases/busco \
                --download prokaryota

            # Info: You can also swap "--download prokaryota" with "--download virus" if you're feeling adventurous ...
            
            touch {output}

            # Clean up 
            rm dummy.fasta
        
        fi

    """



rule checkm2_download:
    output:
        #flag = touch("{base_variable}/databases/checkm2/checkm2_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        database = expand("{base_variable}/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd", base_variable = base_variable),
    params:
        download_path = expand("{base_variable}/databases/checkm2", base_variable = base_variable),
    conda: "conda_definitions/checkm2_conda.yaml"
    shell: """


        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Checking for internet access using google."
            ping -q -c1 google.com &>/dev/null && echo "Online" || echo "Warning: It seems like you don't have internet access? Downloading will probably fail."

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."
        
            checkm2 database \
                --download \
                --path {params.download_path}

            # Consider running checkm2 testrun. Is time and resource consuming though.
            # checkm2 testrun 
            
            touch {output}
        
        fi

    """



rule kraken2_download:
    output:
        flag = touch("{base_variable}/databases/kraken2/hash.k2d") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
    #params: 
        #asscom2_kraken_db = config['asscom2_kraken2_db'],
    conda: "conda_definitions/curl.yaml"
    shell: """


        # TODO: Make a way to check for internet access instead of just crashing. Same for busco and checkm2

        ## Pick a db from this list
        # https://benlangmead.github.io/aws-indexes/k2

        
        ## Shortcuts. Select no bigger than the size of your RAM
    
        #db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz"      # Standard 49GB
        #db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz" # Standard  8GB
        db_pick="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20230314.tar.gz" # Standard 16GB
        

        db_destination="{wildcards.base_variable}/databases/kraken2/kraken2_db.tar.gz"

        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Checking for internet access using google."
            ping -q -c1 google.com &>/dev/null && echo "Online" || echo "Warning: It seems like you don't have internet access? Downloading will probably fail."

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            >&2 echo "Downlading $db_pick to $db_destination"
            mkdir -p $(dirname "$db_destination")
            curl "$db_pick" \
                --output "$db_destination"

            >&2 echo "Decompressing ..."
            tar \
                -xfv $db_destination \
                --directory $(dirname $db_destination)

            >&2 echo "kraken2 DB setup completed"
            echo "Downloaded $db_pick at $(date -Iseconds)" > $(dirname $db_destination)/info.txt
            touch {output}

        fi

    """




rule gtdb_download:
    output:
        #flag = touch("{base_variable}/databases/gtdb/gtdb_download_done.flag") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.
        flag = touch("{base_variable}/databases/gtdb/release207_v2/taxonomy/gtdb_taxonomy.tsv") # Be aware that using snakemake --forcerun will delete the output before rerunning, thus the flag will _always_ be missing. This is  only relevant during development.

    #params: 
        #asscom2_gtdb_db = config['asscom2_gtdb_db'],
    conda: "conda_definitions/curl.yaml"
    shell: """

        # TODO: Make a way to check for internet access instead of just crashing. Same for busco and checkm2

        # https://ecogenomics.github.io/GTDBTk/installing/index.html

        # Pick a source file
        db_pick="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz"
        #db_pick="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz" # alternative mirror, maybe faster in europe? Seems a bit unstable at time of writing.

        db_destination="{wildcards.base_variable}/databases/gtdb/gtdb_db.tar.gz" # Should be defined from 

        # If some previous batch of asscom2 has downloaded the database, we'll just reuse it.
        if [ -f "{output}" ]; then    

            >&2 echo "Flag exists already: touch it to update the mtime ..."
            touch {output}
            
        else

            >&2 echo "Checking for internet access using google."
            ping -q -c1 google.com &>/dev/null && echo "Online" || echo "Warning: It seems like you don't have internet access? Downloading will probably fail."

            >&2 echo "Flag doesn't exist: Download the database and touch the flag ..."

            >&2 echo "Downlading $db_pick to $db_destination"
            mkdir -p $(dirname "$db_destination")
            curl "$db_pick" \
                --output "$db_destination"

            >&2 echo "Decompressing ..."
            tar \
                -xfv $db_destination \
                --directory $(dirname $db_destination)

            >&2 echo "gtdb DB setup completed"
            echo "Downloaded $db_pick at $(date -Iseconds)" > $(dirname $db_destination)/info.txt
            touch {output}

        fi

    """


# ---- Other -------------------------------------------------------------------








# --- Targets for each sample below: --------------------------------

rule sequence_lengths_individual:
    input: "{results_directory}/samples/{sample}/{sample}.fa"
    output: "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    resources:
        runtime = "01:00:00",
        mem_mb = 128,
    conda: "conda_definitions/seqkit.yaml"
    shell: """

        # TODO: Consider whether seqkit stats might be faster?

        seqkit fx2tab {input} -l -g -G -n -H \
        > {output}

    """






rule prokka_individual:
    input: "{results_directory}/samples/{sample}/{sample}.fa"
    output:
        gff = "{results_directory}/samples/{sample}/prokka/{sample}.gff",
        log = "{results_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{results_directory}/samples/{sample}/prokka/{sample}.tsv",
        gff_nofasta = "{results_directory}/samples/{sample}/prokka/{sample}.gff_nofasta",
    container: "docker://staphb/prokka"
    conda: "conda_definitions/prokka.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.prokka_individual.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
      
        # For debugging of minced
        minced --version 
        
        prokka \
            --cpus {threads} \
            --force \
            --rfam \
            --outdir {wildcards.results_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input} \
        | tee {output.log} 

        # I don't remember what I'm actually using this output for?
        # Remove fasta from gff and add sample label
        gff_fasta_start=$(grep --line-number --extended-regexp "^##FASTA" {output.gff} | cut -f1 -d:)
        head --lines $(($gff_fasta_start-1)) {output.gff} \
        > {output.gff_nofasta}

    """











# Kraken is for reads, so why are we using it here without shredding the reads, or at least doing some proportion analysis.


# It looks like loading the database takes as long as analysing  the sample, so it might be faster (and more cpu-hour efficient) to just run all samples in the same job. Unfortunately, it doesn't seem to be possible.
rule kraken2_individual:
    input: 
        assembly = "{results_directory}/samples/{sample}/{sample}.fa",
        database = expand("{base_variable}/databases/kraken2/hash.k2d", base_variable = base_variable),
    output: 
        report = "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv",
        full = "{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_full.tsv",
    params: 
        #asscom2_kraken2_db = config["asscom2_kraken2_db"],
           base_variable = base_variable,
    container: "docker://staphb/kraken2"
    conda: "conda_definitions/kraken2.yaml"
    threads: 2
    resources:
        mem_mb = 16384,
    benchmark: "{results_directory}/benchmarks/benchmark.kraken2_individual.{sample}.tsv"
    shell: """

        db_path="{params.base_variable}/databases/kraken2"
        echo using kraken2 database $db_path


        # Run kraken2
        # https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
        kraken2 \
            --threads {threads} \
            --db $db_path \
            --confidence 0.15 \
            --report {output.report} \
            --report-minimizer-data \
            {input.assembly} \
            > {output.full}

        # Argument on confidence parameter https://www.biostars.org/p/402619/


    """


rule busco_individual:
    input: 
        #busco_download = expand("{base_variable}/databases/busco/busco_download_done.flag", base_variable = base_variable),
        busco_download = expand("{base_variable}/databases/busco/file_versions.tsv", base_variable = base_variable),
        fasta = "{results_directory}/samples/{sample}/{sample}.fa",
    output: 
        flag = touch("{results_directory}/samples/{sample}/busco/busco_done.flag"),
        table_extract = "{results_directory}/samples/{sample}/busco/short_summary_extract.tsv"
    params:
        base_variable = base_variable,
        #results_directory = results_directory,
        out_dir = "{results_directory}/samples/{sample}/busco",
    conda: "conda_definitions/busco.yaml"
    threads: 1 # Because run_sepp hangs for a long time, not doing anything, I'd rather have more processes started on my small CPU.
    resources:
        mem_mb = 8192,
        runtime = "06:00:00",
    shell: """

        # Busco fails because of a problem with the sepp package. This doesn't really matter as we just want the completeness results.
        # But, this means that we need a hacky workaround to let this job exit gracefully (exit code 0) on the basis of whether any completeness results have been written to disk.
        # Hence, the actual exit code of busco, we will ignore.


        >&2 echo "Busco individual"
        # https://busco.ezlab.org/busco_userguide.html#offline
        busco \
            --cpu {threads} \
            --in {input.fasta} \
            --out {params.out_dir} \
            --mode geno \
            --auto-lineage-prok \
            --force \
            --tar \
            --download_path {params.base_variable}/databases/busco \
            --offline || (>&2 echo "ac2: busco failed internally")


        >&2 echo "debug1"
        ### New: Using JSON


        # Cat all auto lineage results together or create empty file
        # The following cat command will fail if the glob doesn't resolve any files: This is the wanted behaviour.
        cat {wildcards.results_directory}/samples/{wildcards.sample}/busco/auto_lineage/*/short_summary.json \
        > {output.table_extract}_temp         
        
        >&2 echo "Results clearly must exist ... "
        
        # Extract relevant features
        cat {output.table_extract}_temp \
        | grep -oE "(\\"in\\"|\\"name\\"|\\"one_line_summary\\").+" \
        > {output.table_extract}

        >&2 echo "debug3"
        # Clean up
        rm {output.table_extract}_temp

        >&2 echo "debug4"

        {void_report}

    """


# --- Collect results among all samples -----------------------------



rule prokka:
    input:
        metadata = expand("{results_directory}/metadata.tsv", results_directory = results_directory),
        gff = expand(
            "{results_directory}/samples/{sample}/prokka/{sample}.gff",
            results_directory = results_directory,
            sample = df["sample"]
        ),


rule busco:
    input: 
        metadata = expand("{results_directory}/metadata.tsv", results_directory = results_directory),
        tables = expand("{results_directory}/samples/{sample}/busco/short_summary_extract.tsv", results_directory = results_directory, sample = df["sample"])


#echo -e "sample\tmatch_percent\tclade_mappings\tlevel_mappings\tlevel\ttaxonomic_id\tclade" \
rule kraken2:
    input: 
        metadata = expand("{results_directory}/metadata.tsv", results_directory = results_directory),
        reports = expand("{results_directory}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", results_directory = results_directory, sample = df["sample"])


rule sequence_lengths:
    input: 
        metadata = expand("{results_directory}/metadata.tsv", results_directory = results_directory),
        lengths = expand("{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", results_directory = results_directory, sample = df["sample"])











# # This one doesn't seem to work, I don't know what is up? Could be nice to have it fixed.
# rule sample_pathway_enrichment_analysis:
#     input: "{results_directory}/collected_results/prokka_labelled.tsv"
#     output: "{results_directory}/collected_results/sample_pathway_enrichment_analysis.tsv"
#     conda: "conda_definitions/r-clusterProfiler.yaml"
#     shell: """


#         Rscript $ASSCOM2_BASE/scripts/sample_pathway_enrichment_analysis.R $ASSCOM2_BASE/assets/ko {input} \
#             > {output}

#         {void_report}
#     """



# --- Targets for the complete set below: ---------------------------

rule checkm2:
    input:
        metadata = "{results_directory}/metadata.tsv",
        #checkm2_download = expand("{base_variable}/databases/checkm2/checkm2_download_done.flag", base_variable = base_variable), #expanding this variable shouldn't be necessary, but it is, because the variable is not present in the output.
        database = expand("{base_variable}/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd", base_variable = base_variable),
        fasta = df["input_file_fasta"].tolist()
    output:
        table = touch("{results_directory}/checkm2/quality_report.tsv"),
    conda: "conda_definitions/checkm2_conda.yaml"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = "24:00:00",
    params:
        rule_dir = results_directory + "/checkm2",
        base_variable = base_variable,
    shell: """

        checkm2 predict \
            --threads {threads} \
            --input {input.fasta} \
            --output-directory {params.rule_dir} \
            --extension .fa \
            --force

        {void_report}

    """



def get_mem_roary(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]

rule roary:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        gff = expand("{results_directory}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], results_directory = results_directory),
    output:
        analyses = ["{results_directory}/roary/summary_statistics.txt", "{results_directory}/roary/core_gene_alignment.aln", "{results_directory}/roary/gene_presence_absence.csv"]
    params:
        blastp_identity = int(config['roary_blastp_identity']), # = 95 # For clustering genes
        core_perc = 99,  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_roary,
        runtime = "23:59:59", # Well, fuck me if this doesn't work on PBS
    container: "docker://sangerpathogens/roary"
    #conda: "conda_definitions/roary.yaml" 
    shell: """
    
        # Since I reinstalled conda, I've had problems with "Can't locate Bio/Roary/CommandLine/Roary.pm in INC". Below is a hacky fix
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0
        

        # Silence parallel's citation pester:
        echo "will cite" | parallel --citation > /dev/null 2> /dev/null

        # Roary is confused by the way snakemake creates directories ahead of time.
        rm -r {wildcards.results_directory}/roary


        roary -a -r -e --mafft \
            -p {threads} \
            -i {params.blastp_identity} \
            -cd {params.core_perc} \
            -f {wildcards.results_directory}/roary \
            {input.gff}
            
        {void_report}

    """


rule snp_dists:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aln = "{results_directory}/roary/core_gene_alignment.aln",
    output: "{results_directory}/snp-dists/snp-dists.tsv"
    conda: "conda_definitions/snp-dists.yaml"
    container: "docker://staphb/snp-dists"
    shell: """

        snp-dists {input.aln} > {output}

        {void_report}
    """




rule assembly_stats:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/assembly-stats/assembly-stats.tsv"
    container: "docker://sangerpathogens/assembly-stats"
    conda: "conda_definitions/assembly-stats.yaml"
    shell: """
        
        assembly-stats -t {input.fasta} > {output}

        {void_report}
    """


def get_mem_gtdbtk(wildcards, attempt): 
    return [150000, 300000, 400000, 500000][attempt-1]

rule gtdbtk:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        db_flag = expand("{base_variable}/databases/gtdb/release207_v2/taxonomy/gtdb_taxonomy.tsv", base_variable = base_variable),
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/gtdbtk/gtdbtk.bac.summary.tsv"
    params:
        batchfile_content = df[['input_file_fasta', 'sample']].to_csv(header = False, index = False, sep = "\t"),
        out_dir = "{results_directory}/gtdbtk/",
        base_variable = base_variable,
        #gtdbtk_data_path = config["gtdbtk_data_path"],
    threads: 8
    #retries: 3
    resources:
        #mem_mb = 150000, # Last time I remember, it used 130000
        mem_mb = get_mem_gtdbtk,
    conda: "conda_definitions/gtdbtk.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gtdbtk.tsv"
    shell: """

        export GTDBTK_DATA_PATH="{params.base_variable}/databases/gtdb/release207_v2" # Should be defined from 

        # Create batchfile
        echo '''{params.batchfile_content}''' > {wildcards.results_directory}/gtdbtk/batchfile.tsv

        gtdbtk classify_wf -h
        
        gtdbtk classify_wf \
            --batchfile {wildcards.results_directory}/gtdbtk/batchfile.tsv \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --keep_intermediates \
            --force

        # Homogenize database version number
        #cp {wildcards.results_directory}/gtdbtk/gtdbtk.bac120.summary.tsv {output}
        # New better version below that also incorporates archaea
        #cat {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv {output}
        

        # Even better: Should be tested on originals
        echo -e "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" \
        > {output}
        tail -n +2 {wildcards.results_directory}/gtdbtk/gtdbtk.*.summary.tsv \
        >> {output}
        


        {void_report}


    """ 



rule abricate:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output:
        card_detailed = "{results_directory}/abricate/card_detailed.tsv",
        card_sum = "{results_directory}/abricate/card_summarized.tsv",
        plasmidfinder_detailed = "{results_directory}/abricate/plasmidfinder_detailed.tsv",
        plasmidfinder_sum = "{results_directory}/abricate/plasmidfinder_summarized.tsv",
        ncbi_detailed = "{results_directory}/abricate/ncbi_detailed.tsv",
        ncbi_sum = "{results_directory}/abricate/ncbi_summarized.tsv",
        vfdb_detailed = "{results_directory}/abricate/vfdb_detailed.tsv",
        vfdb_sum = "{results_directory}/abricate/vfdb_summarized.tsv",

    container: "docker://staphb/abricate"
    conda: "conda_definitions/abricate.yaml"
    shell: """


        # TODO: update these databases

        abricate --db ncbi {input.fasta} > {output.ncbi_detailed}
        abricate --summary {output.ncbi_detailed} > {output.ncbi_sum}

        abricate --db card {input.fasta} > {output.card_detailed}
        abricate --summary {output.card_detailed} > {output.card_sum}
        
        abricate --db plasmidfinder {input.fasta} > {output.plasmidfinder_detailed}
        abricate --summary {output.plasmidfinder_detailed} > {output.plasmidfinder_sum}

        abricate --db vfdb {input.fasta} > {output.vfdb_detailed}
        abricate --summary {output.vfdb_detailed} > {output.vfdb_sum}


        {void_report}
    """



# Parse the mlst scheme for bash
if config["mlst_scheme"] == "automatic":
    mlst_scheme_interpreted = "",
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}",
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{results_directory}/mlst/mlst_schemes.txt",
    container: "docker://staphb/mlst"
    conda: "conda_definitions/mlst.yaml"
    shell: """

        mlst {params.mlst_scheme_interpreted} {input.fasta} > {output}

        mlst --list > {params.list_}

        {void_report}
    """




rule mashtree:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: 
        tree = "{results_directory}/mashtree/mashtree.newick",
        dist = "{results_directory}/mashtree/mash_dist.tsv",
    container: "docker://staphb/mashtree"
    conda: "conda_definitions/mashtree.yaml"
    threads: 4
    resources:
        mem_mb = 16000,
    shell: """

        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist} \
            {input.fasta} > {output.tree}


        {void_report}
    """ 

# TODO:
#rule mash_screen:

def get_mem_fasttree(wildcards, attempt): 
    return [8000, 64000][attempt-1]


rule fasttree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = "{results_directory}/roary/core_gene_alignment.aln",
    output: "{results_directory}/fasttree/fasttree.newick"
    container: "docker://staphb/fasttree"
    conda: "conda_definitions/fasttree.yaml"
    threads: 4
    #retries: 1
    resources:
        mem_mb = get_mem_fasttree,
        runntime = "23:59:59",
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree \
            -nt \
            -gtr {input.fasta} \
        > {output} \
        2> {output}.log 

        {void_report}

    """


rule fetch_report_template:
    output: "{results_directory}/rmarkdown_template.rmd"
    shell: """

        cp $ASSCOM2_BASE/scripts/genomes_to_report_v2.Rmd {output}

        {void_report}
    """



# This rule might seem silly, but it makes sure that the report environment is ready to rock when the report subpipeline eventually is run: This has two pros:
#    1) The vastly faster mamba configuration in the asscom2 pipeline is used
#    2) The conda/mamba debugging is taken care of, without having to wait for jobs to finish on fresh installations.
# Since all snakemake conda environments are installed in $SNAKEMAKE_CONDA_PREFIX set to ${ASSCOM2_BASE}/conda_base, reuse is guaranteed.
rule install_report_environment_aot:
    output: touch(f"{results_directory}/.install_report_environment_aot.flag")
    conda: "report_subpipeline/conda_definitions/r-markdown.yaml"
    shell: """

        echo OK

    """

# Just a dummy rule if you wanna force the report
# assemblycomparator2 --until report
# It isn't enough to just touch the file. The report_subpipeline will not be triggered if the file is empty. Thus we add the date, and we have a nice debug log for seeing when the report was triggered.
# Will only but run if asked to. No need to use --forcerun, since snakemake states this in the output: "reason: Rules with neither input nor output files are always executed."
rule report:
    run: # No need to allocate a job on the cluster for this small job.
        now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        with open(f"{results_directory}/.asscom2_void_report.flag", "a") as void_report_file:
            void_report_file.write(f"ss_{now}\n")



# Call the report subpipeline
report_call = f"""
    mkdir -p {results_directory}/logs; \
    snakemake \
        --snakefile $ASSCOM2_BASE/report_subpipeline/snakefile \
        --cores 4 \
        --use-conda \
        -p \
        --config results_directory=$(pwd)/{results_directory} base_variable={base_variable} batch_title={batch_title} 2> {results_directory}/logs/report.err.log 
    """

onsuccess:
    print("On success: calling report subpipeline ...")
    shell(report_call)

onerror:
    print("On error: calling report subpipeline ...")
    shell(report_call)





print("*/")




