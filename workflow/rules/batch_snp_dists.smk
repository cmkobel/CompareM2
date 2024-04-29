


rule snp_dists:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        #aln = "{results_directory}/roary/core_gene_alignment.aln",
        aln = core_genome_if_exists,
    output: "{results_directory}/snp-dists/snp-dists.tsv"
    conda: "../envs/snp-dists.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.snp_dists.tsv"
    threads: 4
    shell: """

        snp-dists \
            -j {threads} \
            {input.aln:q} > {output:q}

        {void_report}
        
    """
