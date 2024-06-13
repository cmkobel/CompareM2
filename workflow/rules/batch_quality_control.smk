
rule checkm2:
    input:
        metadata = "{output_directory}/metadata.tsv",
        database_representative = DATABASES + "/checkm2/comparem2_checkm2_database_representative.flag",
        fasta = df["input_file_fasta"].tolist()
    output:
        table = touch("{output_directory}/checkm2/quality_report.tsv"),
    conda: "../envs/checkm2.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.checkm2.tsv"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = "24h",
    params:
        rule_dir = output_directory + "/checkm2",
        base_variable = base_variable,
    shell: """

        # Collect version number.
        echo "checkm2 $(checkm2 --version)" > "$(dirname {output})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.table})/.database_version.txt"
        
        checkm2 predict \
            --threads {threads} \
            --input {input.fasta:q} \
            --output-directory {params.rule_dir:q} \
            --database_path $(dirname {input.database_representative:q})/CheckM2_database/uniref100.KO.1.dmnd \
            --force

        {void_report}

    """
