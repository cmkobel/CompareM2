# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent? Or just remove it. I think it creates as many problems as it solves.
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: 
        fasta = "{output_directory}/samples/{sample}/{sample}.fna",
        md5sum = "{output_directory}/samples/{sample}/{sample}.md5.txt",
    #conda: "../envs/any2fasta.yaml"
    threads: 1 # Weirdly, or bugly, there must be a thread n definition in the rule. Otherwise, the set-threads option (in the orion profile) will not be taken up. 
    resources:
        mem_mb = 256,
        runtime = "10m",
    shell: """

        # Collect version number.
        #any2fasta -v > "$(dirname {output.fasta})/.software_version.txt"
        
        #any2fasta {input.genome:q} > {output.fasta:q}
        cp {input.genome:q} {output.fasta:q}
        
        md5sum {output.fasta:q} > {output.md5sum:q}
        
    """  


rule assembly_stats:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{output_directory}/assembly-stats/assembly-stats.tsv"
    conda: "../envs/assembly-stats.yaml"
    benchmark: "{output_directory}/benchmarks/benchmarks.assembly_stats.tsv"
    shell: """
    
        # Collect version number.
        echo "assembly-stats $(assembly-stats -v)" > "$(dirname {output})/.software_version.txt"
        
        assembly-stats \
            -t \
            {input.fasta:q} > {output:q}

        {void_report}

    """


rule sequence_lengths:
    input:
        metadata = "{output_directory}/metadata.tsv",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna", 
    output: "{output_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    threads: 1
    resources:
        runtime = "60m",
    conda: "../envs/seqkit.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.sequence_lengths_sample.{sample}.tsv"
    shell: """
    
        # Collect version number.
        echo "seqkit $(seqkit -h | grep 'Version:')" > "$(dirname {output})/.software_version.txt"
        
        seqkit fx2tab \
            {input.assembly:q} \
            -l \
            -g \
            -G \
            -n \
            -H \
        > {output:q}

    """


