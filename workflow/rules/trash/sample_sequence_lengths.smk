rule sequence_lengths:
    input:
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna", 
    output: "{results_directory}/samples/{sample}/sequence_lengths/{sample}_seqlen.tsv"
    threads: 1
    resources:
        runtime = "60m",
    conda: "../envs/seqkit.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.sequence_lengths_individual.{sample}.tsv"
    shell: """

        seqkit fx2tab {input.assembly:q} -l -g -G -n -H \
        > {output:q}

    """
