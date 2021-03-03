


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
#   any2fasta (wide input format support)
#   prokka (annotation)
#   kraken2 (species identification)
#   mlst (multi locus sequence typing)
#   abricate (virulence/resistance gene identification)
#   (Oriloc) (Identify possible replication origins, and thereby identify chromids)
# For each group
#   roary (pan and core genome)
#   snp-dists (core genome snp-distances)
#   panito (average nucleotide identity
#   FastTree (phylogenetic tree of core genome)
#   IQ-tree (phylogenetic tree of core genome with bootstrapping)
#   (GC3-profiling) ("fingerprinting" of the distribution of GC-content)
#   (Identification of horizontally transferred genes)


title = "test"


out_base_var = "output_asscom1"




#sample_file = "input/" + title + ".tsv"
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


print(df)
print("///")



# Collect all targets
rule all:
    input: expand(["{out_base}/sample_{sample}/{sample}.fa", \
                   "{out_base}/sample_{sample}/prokka/{sample}.gff", \
                   "{out_base}/roary/summary_statistics.txt"], \
                  out_base = out_base_var, sample = df["sample"]) # copy


  




# Copy the input file to its new home
rule copy:
    #input: "{sample}"
    input: lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0]
    output: "{out_base}/sample_{sample}/{sample}.fa"
    #log: "logs/{out_base}_{wildcards.sample}.out.log"

    shell: """

        cp {input} {output}

        """


# Prokka
rule prokka:
    input: "{out_base}/sample_{sample}/{sample}.fa"
    output: "{out_base}/sample_{sample}/prokka/{sample}.gff"
    #conda: "envs/prokka.yml"
    container: "docker://staphb/prokka"
    threads: 4
    shell: """

        prokka --cpus 4 --force --outdir {wildcards.out_base}/sample_{wildcards.sample}/prokka --prefix {wildcards.sample} {input} || echo exit 0

        echo "prokka done"
        """



rule roary:
    input: expand("{out_base}/sample_{sample}/prokka/{sample}.gff", sample = df["sample"], out_base = out_base_var)
    output: "{out_base}/roary/summary_statistics.txt"
    params:
        blastp_identity = 95,
        core_perc = 99
    #conda: "envs/roary.yml"
    container: "docker://sangerpathogens/roary"
    threads: 4
    shell: """

        parallel --citation > /dev/null 2>&1
        roary -a -r -e --mafft -p 4 -f {wildcards.out_base}/roary -i {params.blastp_identity} -cd {params.core_perc} {input}

        # todo: use kraken as well
        

        echo "roary done"
        """






