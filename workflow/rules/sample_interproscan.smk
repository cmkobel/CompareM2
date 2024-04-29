
# By default it runs TIGRFAM, Hamap, Pfam
rule interproscan:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        # No external database is needed.
    output:
        tsv = "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv",
    params:
        file_base = "{results_directory}/samples/{sample}/interproscan/{sample}_interproscan",
    conda: "../envs/interproscan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{results_directory}/benchmarks/benchmark.interproscan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """

        # https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html#command-line-options
        interproscan.sh \
            --applications TIGRFAM,Hamap,Pfam \
            --cpu {threads} \
            --output-file-base {params.file_base:q} \
            --disable-precalc \
            --formats TSV \
            --goterms \
            --iprlookup \
            --pathways \
            --seqtype p \
            --tempdir {resources.tmpdir} \
            --input {input.aminoacid:q}

    """
