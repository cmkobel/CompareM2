  


rule sequence_lengths: # TODO: rename to seqkit?
    input:
        assembly = "{output_directory}/samples/{sample}/{sample}.fna", 
        #metadata = "{output_directory}/metadata.tsv", # Only batch jobs need metadata (for proper updating)
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


