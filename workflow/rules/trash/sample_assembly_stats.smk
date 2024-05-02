rule assembly_stats:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/assembly-stats/assembly-stats.tsv"
    conda: "../envs/assembly-stats.yaml"
    benchmark: "{results_directory}/benchmarks/assembly_stats.tsv"
    shell: """
        
        assembly-stats -t {input.fasta:q} > {output:q}

        {void_report}

    """
