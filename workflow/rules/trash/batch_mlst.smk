
# Parse the mlst scheme for bash
if config["mlst_scheme"] == "automatic":
    mlst_scheme_interpreted = "",
else:
    mlst_scheme_interpreted = f"--scheme {config['mlst_scheme']}",
#print(f"Info: The mlst_scheme is set to <{mlst_scheme_interpreted}>") # Debug message.

rule mlst:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        fasta = df["input_file_fasta"].tolist(),
    output: "{results_directory}/mlst/mlst.tsv",
    params:
        mlst_scheme_interpreted = mlst_scheme_interpreted,
        list_ = "{results_directory}/mlst/mlst_schemes.txt", 
    threads: 4
    conda: "../envs/mlst.yaml"
    benchmark: "{results_directory}/benchmarks/mlst.tsv"
    shell: """

        mlst \
            --threads {threads} {params.mlst_scheme_interpreted} \
            {input.fasta:q} \
            > {output:q}

        # Dump available mlst databases
        mlst --list > {params.list_}

        {void_report}

    """
