
# snakemake --snakefile ~/assemblycomparator2/snakefile --profile ~/assemblycomparator2/configs/slurm/ --cluster-config ~/assemblycomparator2/configs/cluster.yaml 

from os import listdir
from os.path import isfile, join
#import yaml
import pandas as pd
import numpy as np
#import time
#import re
#from shutil import copyfile
#import re


print()
print("         █████╗ ███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗  ")
print("        ██╔══██╗██╔════╝██╔════╝██╔════╝██╔═══██╗████╗ ████║╚════██╗ ")
print("        ███████║███████╗███████╗██║     ██║   ██║██╔████╔██║ █████╔╝ ")
print("        ██╔══██║╚════██║╚════██║██║     ██║   ██║██║╚██╔╝██║██╔═══╝  ")
print("        ██║  ██║███████║███████║╚██████╗╚██████╔╝██║ ╚═╝ ██║███████╗ ")
print("        ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝ ")
print("                      A.K.A. assemblycomparator2                     ")
print()





out_base_var = "output_asscom2"



#reference = config["reference"]


# ---- Read in files in the current working directory ---------------

relative_wd = "."
#extension_whitelist = ["fna", "fa", "fas", "fasta", "seq"] # old 
extension_whitelist = ["fna", "fa", "fas", "fasta", "seq", "gb", "fq", "gff", "gfa", "clw", "sth", "gz", "bz2", "zip"]

present_files = [f for f in listdir(relative_wd) if isfile(join(relative_wd,f))]


df = pd.DataFrame(data = {'input_file': present_files})

# Check that the directory is not empty.
if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...")
    raise Exception("Zero files.")



df = df[~df["input_file"].str.startswith(".", na = False)] # Remove hidden files
df['sample'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()] # Extract everything before the extension dot.
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()] # Extract extension
df['input_file_fasta'] = out_base_var + "/samples/" + df['sample'] + "/" + df['sample'] + ".fa" # This is where the input file is copied to in the first snakemake rule.



  
# ---- Displaying filtered dataframe ready for analysis -------------

df = df.reset_index(drop = True)
print(df)
print("//")
print()


# Ask the user if they want to continue
#input_continue = input("continue? (y/n) ")
#if not input_continue.lower()[0] == "y":
#    exit("Quitting ...")
#


# Create the dirs

try:
    os.mkdir(out_base_var)
    os.mkdir("logs") # The log directory is actually not used for local setups
except OSError as e:
    print(e)
else:
    print ("Successfully created the output directories.")


# Collect all targets
rule all:
    input: expand(["{out_base}/metadata.tsv", \
                   "{out_base}/samples/{sample}/{sample}.fa", \
                   "{out_base}/assembly-stats/assembly-stats.tsv", \
                   "{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", \
                   "{out_base}/samples/{sample}/prokka/{sample}.gff", \
                   "{out_base}/kraken2/kraken2_reports.tsv", \
                   "{out_base}/roary/summary_statistics.txt", \
                   "{out_base}/abricate/card_detailed.tsv", \
                   "{out_base}/mlst/mlst.tsv", \
                   "{out_base}/mashtree/mashtree.newick", \
                   "{out_base}/fasttree/fasttree.newick", \
                   "{out_base}/report.html"], \
                  out_base = out_base_var, sample = df["sample"]) # copy


  

# Write the df table to the directory for later reference.
rule metadata:
    input: df["input_file"].tolist()
    output: "{out_base}/metadata.tsv"
    run: 
        df.to_csv(str(output), index_label = "index", sep = "\t")



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

        #cp {input} {output}

        any2fasta {input} > {output}


    """

# rule assembly_stats:
#     input: "{out_base}/samples/{sample}/{sample}.fa"
#     output: "{out_base}/samples/{sample}/assembly-stats/{sample}_assemblystats.txt"
#     conda: "conda_envs/assembly-stats.yaml"
#     shell: """
        
#         assembly-stats -t {input} > {output}
    
#     """




##################################
# Targets for each sample below: #
##################################


rule seqlen:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    container: "docker://staphb/bioawk"
    conda: "conda_envs/bioawk.yaml"
    shell: """

        bioawk -v sam={wildcards.sample} -c fastx '{{ print sam, $name, length($seq) }}' < {input} \
        > {output}

    """



rule prokka:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output:
        gff = "{out_base}/samples/{sample}/prokka/{sample}.gff",
        txt = "{out_base}/samples/{sample}/prokka/{sample}.txt",
        summary = "{out_base}/samples/{sample}/prokka/{sample}_summary.txt",


    #conda: "envs/prokka.yml"
    container: "docker://staphb/prokka"
    conda: "conda_envs/prokka.yaml"
    threads: 4
    shell: """

        prokka --cpus {threads} --force --outdir {wildcards.out_base}/samples/{wildcards.sample}/prokka --prefix {wildcards.sample} {input} #|| echo exit 0

        # Put sample names in front
        cat {output.txt} \
        | awk -v sam={wildcards.sample} '{{ print sam ": " $0 }}' \
        > {output.summary}



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


# rule parse_kraken2:
#     input: "{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv"
#     output: "{out_base}/samples/{sample}/kraken2/{sample}_kraken2_top10.tsv"
#     run:
#         print("hej")
#         kraken2_report = pd.read_csv(str(input),
#             sep = '\t',
#             names = ['match_percent', 'clade_mappings', 'level_mappings', 'level', 'taxonomic_id', 'clade'],
#             dtype = str)

#         kraken2_report.insert(loc=0, column='sample', value=wildcards.sample)

#         # Remove superfluous spaces
#         kraken2_report = kraken2_report.apply(np.vectorize(lambda x: str(x).strip()))

#         print(kraken2_report)


#         kraken2_report.to_csv(str(output), index = False, sep = "\t")



rule collect_kraken2:
    input:
        kraken2 = expand("{out_base}/samples/{sample}/kraken2/{sample}_kraken2_report.tsv", out_base = out_base_var, sample = df["sample"]),
        seqlen = expand("{out_base}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv", out_base = out_base_var, sample = df["sample"])

    output:
        kraken2 = "{out_base}/kraken2/kraken2_reports.tsv",
        seqlen = "{out_base}/sequence_lengths/sequence_lengths.tsv"

    shell: """

        # kraken2
        echo -e "sample\tmatch_percent\tclade_mappings\tlevel_mappings\tlevel\ttaxonomic_id\tclade" \
        > {output.kraken2}

        cat {input.kraken2} >> {output.kraken2}


        # Sequence lengths
        cat {input.seqlen} > {output.seqlen} 



    """





#######################################
# Targets for the complete set below: #
#######################################
rule roary:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: ["{out_base}/roary/summary_statistics.txt", "{out_base}/roary/core_gene_alignment.aln", "{out_base}/roary/gene_presence_absence.csv", "{out_base}/roary/roary_done.flag"]
    params:
        blastp_identity = 95, # For clustering genes
        core_perc = 99  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 4
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



rule assembly_stats:
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/assembly-stats/assembly-stats.tsv"
    conda: "conda_envs/assembly-stats.yaml"
    shell: """
        
        assembly-stats -t {input} > {output}
    
    """




rule abricate:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    #input: df["input_file"].tolist()
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


rule mlst:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/mlst/mlst.tsv"
    container: "docker://staphb/mlst"
    conda: "conda_envs/mlst.yaml"
    shell: """

        mlst {input} > {output}

    """




rule mashtree:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file_fasta"].tolist()
    output: "{out_base}/mashtree/mashtree.newick"
    container: "docker://staphb/mashtree"
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

        FastTree -nt -gtr {input} > {output} 2> {output}.log

    """



# rule roary_plots:
#     input: genes = "{out_base}/roary/gene_presence_absence.csv",
#         tree = "{out_base}/fasttree/fasttree.newick"
#     output: "{out_base}/roary_plots/whatever"
#     container: "docker://python" # Make our own python container with cairosvg perl etc...
#     shell: """
        
#         # Failing because matplotlib is missing...
#         python3 scripts/roary_plots.py {input.tree} {input.genes} > hat 2> hat.err

#         # TODO: add the other weird stuff from https://github.com/cmkobel/assemblycomparator/blob/61c9a891a75e2f252dc54185d74c0fbb092815e5/workflow_templates.py#L489
#     """






rule report:
    input: "{out_base}/roary/roary_done.flag"
    output: "{out_base}/report.html"
    params:
        markdown_template_rmd = "genomes_to_report_v2.Rmd",
        markdown_template_html = "genomes_to_report_v2.html"
    singularity: "docker://marcmtk/sarscov2_markdown"
    conda: "conda_envs/r-markdown.yaml"
    shell: """

        cd {wildcards.out_base}

        cp $ASSCOM2_BASE/scripts/{params.markdown_template_rmd} .

        Rscript -e 'library(rmarkdown); paste("t", getwd()); rmarkdown::render("{params.markdown_template_rmd}", "html_document")'

        rm {params.markdown_template_rmd}
        mv {params.markdown_template_html} ../{output}
        
    """


#print(mashtree.input)







