
rule iqtree:
    input:
        metadata = "{results_directory}/metadata.tsv",
        fasta = core_genome_if_exists,
    output: 
        newick = "{results_directory}/iqtree/core_genome_iqtree.treefile"
    params:
        version_file = "{results_directory}/iqtree/iqtree_ac2_version.txt"
    conda: "../envs/iqtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.iqtree.tsv"
    threads: 16
    retries: 3
    resources:
        mem_mb = 32000,
        runtime = "24h",
    shell: """

        iqtree --version > {params.version_file}

        iqtree \
            -s {input.fasta:q} \
            -m GTR \
            --boot 100 \
            --prefix $(dirname {output.newick:q})/core_genome_iqtree \
            -redo

        # {void_report} Not in the report yet.


    """
