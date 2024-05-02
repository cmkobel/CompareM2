
def get_mem_fasttree(wildcards, attempt): 
    return [16000, 32000, 64000, 0][attempt-1]
rule fasttree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: "{results_directory}/fasttree/fasttree.newick"
    conda: "../envs/fasttree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.fasttree.tsv"
    threads: 4
    retries: 2
    resources:
        mem_mb = get_mem_fasttree,
        runtime = "24h",
    shell: """

        OMP_NUM_THREADS={threads}

        FastTree \
            -nt \
            -gtr {input.fasta:q} \
        > {output:q} \
        2> {output:q}.log 

        {void_report}

    """
