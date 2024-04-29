# Copy the input file to its new home
# Homogenizes the file extension as well (.fa)
# Should I rename this to "rule any2fasta" just to make it more transparent? 
rule copy:
    input: 
        genome = lambda wildcards: df[df["sample"]==wildcards.sample]["input_file"].values[0],
    output: "{results_directory}/samples/{sample}/{sample}.fna"
    conda: "../envs/any2fasta.yaml"
    threads: 1 # Weirdly, or bugly, there must be a thread n definition in the rule. Otherwise, the set-threads option (in the orion profile) will not be taken up. 
    resources:
        mem_mb = 256,
        runtime = "10m",
    shell: """


        
        any2fasta {input.genome:q} > {output:q}

    """  
