rule mashtree:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: 
        tree = "{results_directory}/mashtree/mashtree.newick",
        dist = "{results_directory}/mashtree/mash_dist.tsv",
    conda: "../envs/mashtree.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.mashtree.tsv"
    threads: 16
    resources:
        mem_mb = 16000,
    shell: """
    
        mashtree \
            --numcpus {threads} \
            --outmatrix {output.dist:q} \
            {input.fasta:q} > {output.tree:q}

        {void_report}
    """ 
