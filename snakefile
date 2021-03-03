


from os import listdir
from os.path import isfile, join
import yaml
import pandas as pd
#import time
#import re
#from shutil import copyfile
#import re




# Roadmap:
#
# For each assembly
#   PEND    any2fasta (wide input format support)
#   OK      prokka (annotation)
#   PEND    kraken2 (species identification)
#   OK      mlst (multi locus sequence typing)
#   OK      abricate (virulence/resistance gene identification)
#   PEND    (Oriloc) (Identify possible replication origins, and thereby identify chromids)
# For each group
#   INP     roary (pan and core genome)
#   PEND    snp-dists (core genome snp-distances)
#   PEND    panito (average nucleotide identity
#   PEND    FastTree (phylogenetic tree of core genome)
#   PEND    IQ-tree (phylogenetic tree of core genome with bootstrapping)
#   PEND    (GC3-profiling) ("fingerprinting" of the distribution of GC-content)
#   PEND    (Identification of horizontally transferred genes)


title = "test"


out_base_var = "output_asscom1"




#samples/file = "input/" + title + ".tsv"
#reference = config["reference"]


# First comes the machinery that lists the assembly files.
mypath = "."
extension_whitelist = ["fna", "fa", "fas", "fasta", "seq"]
present_files = [f for f in listdir(mypath) if isfile(join(mypath,f))]


df = pd.DataFrame(data = {'input_file': present_files})

df['sample'] = [".".join(i.split(".")[:-1]) for i in df['input_file'].tolist()]
df['extension'] =  [i.split(".")[-1] for i in df['input_file'].tolist()]


df = df.loc[df['extension'].isin(extension_whitelist)]

if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...")
    exit("Zero files.")

df = df.reset_index(drop = True)
print(df)
print("//")


#input_continue = input("continue? (y/n) ")
#if not input_continue.lower()[0] == "y":
#    exit("Quitting ...")
#


# Collect all targets
rule all:
    input: expand(["{out_base}/metadata.tsv", \
                   "{out_base}/samples/{sample}/{sample}.fa", \
                   "{out_base}/samples/{sample}/prokka/{sample}.gff", \
                   "{out_base}/roary/summary_statistics.txt", \
                   "{out_base}/abricate/card_detail.tsv", \
                   "{out_base}/mlst/mlst.tsv"], \
                  out_base = out_base_var, sample = df["sample"]) # copy


  

# Write the df table to the directory for later reference.
rule metadata:
    input: df["input_file"].tolist()
    output: "{out_base}/metadata.tsv"
    run: 
        df.to_csv(str(output), index_label = "index", sep = "\t")


# Copy the input file to its new home
rule copy:
    #input: "{sample}"
    input: lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0]
    output: "{out_base}/samples/{sample}/{sample}.fa"
    #log: "logs/{out_base}_{wildcards.sample}.out.log"

    shell: """

        cp {input} {output}

        """


##################################
# Targets for each sample below: #
##################################
rule prokka:
    input: "{out_base}/samples/{sample}/{sample}.fa"
    output: "{out_base}/samples/{sample}/prokka/{sample}.gff"
    #conda: "envs/prokka.yml"
    container: "docker://staphb/prokka"
    threads: 4
    shell: """

        prokka --cpus {threads} --force --outdir {wildcards.out_base}/samples/{wildcards.sample}/prokka --prefix {wildcards.sample} {input} || echo exit 0

        echo "prokka done"
        """


#######################################
# Targets for the complete set below: #
#######################################
rule roary:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: "{out_base}/roary/summary_statistics.txt"
    params:
        blastp_identity = 95,
        core_perc = 99
    #conda: "envs/roary.yml"
    threads: 8
    container: "docker://sangerpathogens/roary"
    shell: """

        roary -a -r -e --mafft -p {threads} -f {wildcards.out_base}/roary -i {params.blastp_identity} -cd {params.core_perc} {input}

        # todo: use kraken as well
        

        """


rule abricate:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output:
        card_detail = "{out_base}/abricate/card_detail.tsv",
        card_sum = "{out_base}/abricate/card_summary.tsv",
        plasmidfinder_detail = "{out_base}/abricate/plasmidfinder_detail.tsv",
        plasmidfinder_sum = "{out_base}/abricate/plasmidfinder_summary.tsv",
        ncbi_detail = "{out_base}/abricate/ncbi_detail.tsv",
        ncbi_sum = "{out_base}/abricate/ncbi_summary.tsv"
    container: "docker://staphb/abricate"
    shell: """


        abricate --db card {input} > {output.card_detail}
        abricate --summary {output.card_detail} > {output.card_sum}
        
        abricate --db plasmidfinder {input} > {output.plasmidfinder_detail}
        abricate --summary {output.plasmidfinder_detail} > {output.plasmidfinder_sum}
        
        abricate --db ncbi {input} > {output.ncbi_detail}
        abricate --summary {output.ncbi_detail} > {output.ncbi_sum}
        


        """




rule mlst:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output: "{out_base}/mlst/mlst.tsv"
    container: "docker://staphb/mlst"
    shell: """

        mlst {input} > {output}

        """






