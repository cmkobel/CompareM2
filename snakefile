
# snakemake --snakefile ~/assemblycomparator2/snakefile --profile ~/assemblycomparator2/configs/slurm/ --cluster-config ~/assemblycomparator2/configs/cluster.yaml 

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
#   PEND    Mashtree
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
# TODO: Filter out hidden files (starting with dot .)

if df.shape[0] == 0:
    print("Error: No fasta files in the current directory. Quitting ...")
    exit("Zero files.")

df = df.reset_index(drop = True)
print(df)
print("//")
print()

#input_continue = input("continue? (y/n) ")
#if not input_continue.lower()[0] == "y":
#    exit("Quitting ...")
#


# Create the dirs

try:
    os.mkdir(out_base_var)
    os.mkdir("logs")
except OSError:
    print ("Creation of the directories")
else:
    print ("Successfully created the directories")


# Collect all targets
rule all:
    input: expand(["{out_base}/metadata.tsv", \
                   "{out_base}/samples/{sample}/{sample}.fa", \
                   "{out_base}/samples/{sample}/prokka/{sample}.gff", \
                   "{out_base}/roary/summary_statistics.txt", \
                   "{out_base}/abricate/card_detail.tsv", \
                   "{out_base}/mlst/mlst.tsv", \
                   "{out_base}/mashtree/mashtree.newick", \
                   "{out_base}/fasttree/fasttree.newick"], \
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

        mkdir -p logs output_asscom1

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

        """


#######################################
# Targets for the complete set below: #
#######################################
rule roary:
    input: expand("{out_base}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: ["{out_base}/roary/summary_statistics.txt", "{out_base}/roary/core_gene_alignment.aln", "{out_base}/roary/gene_presence_absence.csv"]
    params:
        blastp_identity = 95, # For clustering genes
        core_perc = 99  # Definition of the core genome
    #conda: "envs/roary.yml"
    threads: 8
    container: "docker://sangerpathogens/roary"
    shell: """


        # Roary is confused by the way snakemake creates directories ahead of time.
        # So I will delete it manually here before calling roary.
        rm -r {wildcards.out_base}/roary

        roary -a -r -e --mafft -p {threads} -i {params.blastp_identity} -cd {params.core_perc} -f {wildcards.out_base}/roary {input}
                
        
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


        # TODO: update these databases

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




rule mashtree:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: df["input_file"].tolist()
    output: "{out_base}/mashtree/mashtree.newick"
    container: "docker://staphb/mashtree"
    threads: 4
    shell: """

        mashtree --numcpus {threads} {input} > {output}

        """




rule fasttree:
    #input: expand("{out_base}/samples/{sample}/{sample}.fa", sample = df["sample"], out_base = out_base_var)
    input: "{out_base}/roary/core_gene_alignment.aln"
    output: "{out_base}/fasttree/fasttree.newick"
    container: "docker://staphb/fasttree"
    threads: 4
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree -nt {input} > {output} 2> {output}.log

        """



rule roary_plots:
    input: genes = "{out_base}/roary/gene_presence_absence.csv",
        tree = "{out_base}/fasttree/fasttree.newick"
    output: "{out_base}/roary_plots/whatever"
    container: "docker://python" # Make our own python container with cairosvg perl etc...
    shell: """
        
        # Failing because matplotlib is missing...
        python3 scripts/roary_plots.py {input.tree} {input.genes} > hat 2> hat.err

        # TODO: add the other weird stuff from https://github.com/cmkobel/assemblycomparator/blob/61c9a891a75e2f252dc54185d74c0fbb092815e5/workflow_templates.py#L489
        """



#print(mashtree.input)







