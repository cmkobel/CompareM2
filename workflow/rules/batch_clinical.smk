
rule abricate:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output:
        ncbi_detailed = "{output_directory}/abricate/ncbi_detailed.tsv",
        ncbi_sum = "{output_directory}/abricate/ncbi_summarized.tsv",
        card_detailed = "{output_directory}/abricate/card_detailed.tsv",
        card_sum = "{output_directory}/abricate/card_summarized.tsv",
        plasmidfinder_detailed = "{output_directory}/abricate/plasmidfinder_detailed.tsv",
        plasmidfinder_sum = "{output_directory}/abricate/plasmidfinder_summarized.tsv",
        vfdb_detailed = "{output_directory}/abricate/vfdb_detailed.tsv",
        vfdb_sum = "{output_directory}/abricate/vfdb_summarized.tsv",
    conda: "../envs/abricate.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.abricate.tsv"
    shell: """

        # Collect version number.
        abricate -v > "$(dirname output.ncbi_detailed)/.software_version.txt"

        abricate --db ncbi {input.fasta:q} > {output.ncbi_detailed:q}
        abricate --summary {output.ncbi_detailed:q} > {output.ncbi_sum:q}

        abricate --db card {input.fasta:q} > {output.card_detailed:q}
        abricate --summary {output.card_detailed:q} > {output.card_sum:q}
        
        abricate --db plasmidfinder {input.fasta:q} > {output.plasmidfinder_detailed:q}
        abricate --summary {output.plasmidfinder_detailed:q} > {output.plasmidfinder_sum:q}

        abricate --db vfdb {input.fasta:q} > {output.vfdb_detailed:q}
        abricate --summary {output.vfdb_detailed:q} > {output.vfdb_sum:q}

        {void_report}
    """


# TODO extract function
# Parse the mlst scheme for bash
if config["mlst_scheme"] == "automatic":
    mlst_scheme_interpreted = "",
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}",
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{output_directory}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{output_directory}/mlst/mlst_schemes.txt", 
    threads: 4
    conda: "../envs/mlst.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.mlst.tsv"
    shell: """
    
        # Collect version number.
        mlst -v > "$(dirname output.ncbi_detailed)/.software_version.txt"

        mlst \
            --threads {threads} {params.mlst_scheme_interpreted} \
            {input.fasta:q} \
            > {output:q}

        # Dump available mlst databases
        mlst --list > {params.list_}

        {void_report}

    """
