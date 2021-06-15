
# snakemake --snakefile ~/assemblycomparator2/snakefile --profile ~/assemblycomparator2/configs/slurm/ --cluster-config ~/assemblycomparator2/configs/cluster.yaml 

__version__ = "v2.1.0"
__author__ = ['Carl M. Kobel', 'Oliver Kjærlund Hansen']

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
cwd = os.getcwd()
batch_title = cwd.split("/")[-1]
print("/*")
print()
print("         █████╗ ███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗  ")
print("        ██╔══██╗██╔════╝██╔════╝██╔════╝██╔═══██╗████╗ ████║╚════██╗ ")
print("        ███████║███████╗███████╗██║     ██║   ██║██╔████╔██║ █████╔╝ ")
print("        ██╔══██║╚════██║╚════██║██║     ██║   ██║██║╚██╔╝██║██╔═══╝  ")
print("        ██║  ██║███████║███████║╚██████╗╚██████╔╝██║ ╚═╝ ██║███████╗ ")
print("        ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝ ")
print("                      A.K.A. assemblycomparator2                     ")
print()
print(f"            batch_title:            {batch_title}")
print(f"            roary_blastp_identity:  {config['roary_blastp_identity']} (default 95)")
print(f"            mlst_scheme:            {config['mlst_scheme']} (default automatic)")
print()



out_base_var = "output_asscom2"


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
    raise Exception("Zero files.")



df = df[~df["input_file"].str.startswith(".", na = False)] # Remove hidden files
df['sample_raw'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()] # Extract everything before the extension dot.
df['sample'] = df['sample_raw'].str.replace(' ','_')
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()] # Extract extension
df['input_file_fasta'] = out_base_var + "/samples/" + df['sample'] + "/" + df['sample'] + ".fa" # This is where the input file is copied to in the first snakemake rule.

df = df[df['extension'].isin(extension_whitelist)] # Remove files with unsupported formats.

# Check that the directory is not empty, again.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...(2)")
    raise Exception("Zero files.")


#df_mini = df_mini.apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # strip whitespace and replace spaces with underscores.

  
# --- Displaying filtered dataframe ready for analysis --------------

df = df.reset_index(drop = True)
print(df)
print("//")
print()




# --- Make sure the log directory exists. ---------------------------
try:
    os.mkdir("logs") # The log directory is actually not used for local setups
except:
    pass

try: 
    os.mkdir("output_asscom2")
except:
    pass






# --- Collect all targets. ------------------------------------------
rule all:
    input: expand(["{out_base}/metadata.tsv", \
                   "{out_base}/assembly-stats/assembly-stats.tsv", \
                   "{out_base}/collected_results/sequence_lengths.tsv", \
                   "{out_base}/collected_results/prokka_summaries.txt", \
                   "{out_base}/collected_results/kraken2_reports.tsv", \
                   "{out_base}/roary/summary_statistics.txt", \
                   "{out_base}/abricate/card_detailed.tsv", \
                   "{out_base}/mashtree/mashtree.newick", \
                   "{out_base}/mlst/mlst.tsv", \
                   "{out_base}/fasttree/fasttree.newick", \
                   "{out_base}/report_{batch_title}.html", \
                   "{out_base}/snp-dists/snp-dists.tsv"], \
                  out_base = out_base_var, sample = df["sample"], batch_title = batch_title) # copy


  

# Write the df table to the directory for later reference.
rule metadata:
    input: expand("{out_base}/samples/{sample}/{sample}.fa", out_base = out_base_var, sample = df["sample"])
    output: "{out_base}/metadata.tsv"
    run: 
        df.to_csv(str(output), index_label = "index", sep = "\t")
        #os.system(f"cp ${{ASSCOM2_BASE}}/scripts/{report_template_file_basename} {out_base_var}")



# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
rule copy:
    #input: "{sample}"
    input: lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0]
    output: "{out_base}/samples/{sample}/{sample}.fa"
    #log: "logs/{out_base}_{wildcards.sample}.out.log"
    container: "docker://pvstodghill/any2fasta"
    conda: "conda_envs/any2fasta.yaml"
    shell: """

        mkdir -p logs {wildcards.out_base}

        ls -l "{input}"

        any2fasta "{input}" > {output}


    """



# rule assembly_stats:
#     input: "{out_base}/samples/{sample}/{sample}.fa"
#     output: "{out_base}/samples/{sample}/assembly-stats/{sample}_assemblystats.txt"
#     conda: "conda_envs/assembly-stats.yaml"
#     shell: """
        
#         assembly-stats -t {input} > {output}
    
#     """





# --- Targets for each sample below: --------------------------------

rule seqlen:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    container: "docker://cmkobel/bioawk"
    conda: "conda_envs/bioawk.yaml"
    shell: """

        bioawk -v sam={wildcards.sample} -c fastx '{{ print sam, $name, length($seq) }}' < {input} \
        > {output}

    """



rule prokka:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output:
        gff = "{out_base}/samples/{sample}/prokka/{sample}.gff",
        log = "{out_base}/samples/{sample}/prokka/{sample}.log",
        summary = "{out_base}/samples/{sample}/prokka/{sample}_summary.txt"
    container: "docker://staphb/prokka"
    conda: "conda_envs/prokka.yaml"
    threads: 4
    shell: """

        prokka \
            --cpus {threads} \
            --force \
            --outdir {wildcards.out_base}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input} \
            > {output.log} #|| echo exit 0

        cat {output.log} \
            | grep "Found" \
            | grep -E "tRNAs|rRNAs|CRISPRs|CDS|unique" \
            | cut -d" " -f 3,4 \
            | awk -v sam={wildcards.sample} '{{ print sam " " $0 }}' \
            >> {output.summary}

    """



rule kraken2:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv"
    container: "docker://staphb/kraken2"
    conda: "conda_envs/kraken2.yaml"
    threads: 4
    shell: """


        if [ ! -z $ASSCOM2_KRAKEN2_DB ]; then
            echo using kraken2 database $ASSCOM2_KRAKEN2_DB

            # Run kraken2
            kraken2 \
                --threads 4 \
                --db $ASSCOM2_KRAKEN2_DB \
                --report {output}_tmp \
                {input} \
                > /dev/null

            # Put sample names in front
            cat {output}_tmp \
            | awk -v sam={wildcards.sample} '{{ print sam "\t" $0 }}' \
            > {output}

            # Remove temp file
            rm {output}_tmp

        else
            echo "The ASSCOM2_KRAKEN2_DB variable is not set, and thus the kraken2 rule and its jobs will not be run. Consider using the scripts/set_up_kraken2.sh script for downloading and linking the latest kraken2 database."
        fi

    """




# --- Collect results among all samples -----------------------------

rule collect_kraken2:
    input: expand("{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", out_base = out_base_var, sample = df["sample"]),
    output: "{out_base}/collected_results/kraken2_reports.tsv",
    shell: """

        # kraken2
        echo -e "sample\tmatch_percent\tclade_mappings\tlevel_mappings\tlevel\ttaxonomic_id\tclade" \
        > {output}

        cat {input} >> {output} 

    """


rule collect_seqlen:
    input: expand("{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", out_base = out_base_var, sample = df["sample"])
    output: "{out_base}/collected_results/sequence_lengths.tsv"
    shell: """

        # Sequence lengths
        echo -e "sample\trecord\tlength" \
        > {output}

        cat {input} >> {output} 

    """


rule collect_prokka:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}_summary.txt", out_base = out_base_var, sample = df["sample"]),
    output: "{out_base}/collected_results/prokka_summaries.txt",
    shell: """

        # prokka
        echo "sample value name" \
        > {output}

        cat {input} >> {output}

    """









# --- Targets for the complete set below: ---------------------------
rule roary:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: ["{out_base}/roary/summary_statistics.txt", "{out_base}/roary/core_gene_alignment.aln", "{out_base}/roary/gene_presence_absence.csv", "{out_base}/roary/roary_done.flag"]
    params:
        blastp_identity = int(config['roary_blastp_identity']), # = 95 # For clustering genes
        core_perc = 99  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 16
    container: "docker://sangerpathogens/roary"
    conda: "conda_envs/roary.yaml"
    shell: """


        # Roary is confused by the way snakemake creates directories ahead of time.
        # So I will delete it manually here before calling roary.
        rm -r {wildcards.out_base}/roary

        roary -a -r -e --mafft \
            -p {threads} \
            -i {params.blastp_identity} \
            -cd {params.core_perc} \
            -f {wildcards.out_base}/roary \
            {input} || echo roary failed

        touch {output}
                
        
    """


rule snp_dists:
    input: "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/snp-dists/snp-dists.tsv"
    conda: "conda_envs/snp-dists.yaml"
    container: "docker://staphb/snp-dists"
    shell: """

        snp-dists {input} > {output}

    """




rule assembly_stats:
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/assembly-stats/assembly-stats.tsv"
    container: "docker://sangerpathogens/assembly-stats"
    conda: "conda_envs/assembly-stats.yaml"
    shell: """
        
        assembly-stats -t {input} > {output}
    
    """




rule abricate:
    input: df["input_file_fasta"].tolist()
    output:
        card_detailed = "{out_base}/abricate/card_detailed.tsv",
        card_sum = "{out_base}/abricate/card_summarized.tsv",
        plasmidfinder_detailed = "{out_base}/abricate/plasmidfinder_detailed.tsv",
        plasmidfinder_sum = "{out_base}/abricate/plasmidfinder_summarized.tsv",
        ncbi_detailed = "{out_base}/abricate/ncbi_detailed.tsv",
        ncbi_sum = "{out_base}/abricate/ncbi_summarized.tsv",
        vfdb_detailed = "{out_base}/abricate/vfdb_detailed.tsv",
        vfdb_sum = "{out_base}/abricate/vfdb_summarized.tsv"

    container: "docker://staphb/abricate"
    conda: "conda_envs/abricate.yaml"
    shell: """


        # TODO: update these databases

        abricate --db card {input} > {output.card_detailed}
        abricate --summary {output.card_detailed} > {output.card_sum}
        
        abricate --db plasmidfinder {input} > {output.plasmidfinder_detailed}
        abricate --summary {output.plasmidfinder_detailed} > {output.plasmidfinder_sum}
        
        abricate --db ncbi {input} > {output.ncbi_detailed}
        abricate --summary {output.ncbi_detailed} > {output.ncbi_sum}

        abricate --db vfdb {input} > {output.vfdb_detailed}
        abricate --summary {output.vfdb_detailed} > {output.vfdb_sum}


    """



# Parse the mlst scheme for bash
if config["mlst_scheme"] == "automatic":
    mlst_scheme_interpreted = ""
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}"
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{out_base}/mlst/mlst_schemes.txt"
    container: "docker://staphb/mlst"
    conda: "conda_envs/mlst.yaml"
    shell: """

        mlst {params.mlst_scheme_interpreted} {input} > {output}

        mlst --list > {params.list_}



    """



rule mashtree:
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/mashtree/mashtree.newick"
    container: "docker://staphb/mashtree"
    conda: "conda_envs/mashtree.yaml"
    threads: 4
    shell: """

        mashtree --numcpus {threads} {input} > {output}

    """ 

# TODO:
#rule mash_screen:



rule fasttree:
    input: "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/fasttree/fasttree.newick"
    container: "docker://staphb/fasttree"
    conda: "conda_envs/fasttree.yaml"
    threads: 4
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree -nt -gtr {input} > {output} 2> {output}.log || echo "fasttree failed"

        touch {output}

    """


rule fetch_report_template:
    output: "{out_base}/rmarkdown_template.rmd"
    shell: """

        cp $ASSCOM2_BASE/scripts/genomes_to_report_v2.Rmd {output}

    """


rule report:
    input:
        roary = "{out_base}/roary/roary_done.flag", # fasttree depends on roary, so the roary dependency is not necessary.
        fasttree = "{out_base}/fasttree/fasttree.newick", 
        snp_dists = "{out_base}/snp-dists/snp-dists.tsv",
        rmarkdown_template = "{out_base}/rmarkdown_template.rmd"
    #output: "{out_base}/report.html"
    output: "{out_base}/report_{batch_title}.html"
    params:
        #markdown_template_rmd = "rmarkdown_template.rmd", # "genomes_to_report_v2.Rmd"
        markdown_template_html = "genomes_to_report_v2.html"
    container: "docker://cmkobel/assemblycomparator2_report"
    conda: "conda_envs/r-markdown.yaml"
    shell: """

        cd {wildcards.out_base}

        Rscript -e 'library(rmarkdown); rmarkdown::render("rmarkdown_template.rmd", "html_document")'

        rm rmarkdown_template.rmd
        mv rmarkdown_template.html ../{output}
        
    """

print("*/")







