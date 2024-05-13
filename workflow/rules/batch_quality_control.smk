
rule checkm2:
    input:
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
        fasta = df["input_file_fasta"].tolist()
    output:
        table = touch("{results_directory}/checkm2/quality_report.tsv"),
    conda: "../envs/checkm2.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.checkm2.tsv"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = "24h",
    params:
        rule_dir = results_directory + "/checkm2",
        base_variable = base_variable,
    shell: """

        # Collect version number.
        echo "checkm2 $(checkm2 --version)" > "$(dirname {output})/.software_version.txt"
        
        checkm2 predict \
            --threads {threads} \
            --input {input.fasta:q} \
            --output-directory {params.rule_dir:q} \
            --extension .fa \
            --database_path $(dirname {input.database_representative:q})/CheckM2_database/uniref100.KO.1.dmnd \
            --force

        {void_report}

    """
