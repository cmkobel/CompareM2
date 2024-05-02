
# Use the same database as checkm2, but run on amino-acid files instead of dna.
# This will be used for a subsequent pathway enrichment analysis.
# Idea for speed up: concatenate all genomes together first, like they do in checkm2. Then we only need to load the database once.
rule diamond_kegg: # or uniref_ko?
    input: 
        metadata = "{results_directory}/metadata.tsv", # For the report
        aminoacid = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # From prokka
        #database_representative = base_variable + "/databases/checkm2/ac2_checkm2_database_representative.flag",
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
        
    output:
        tsv = "{results_directory}/samples/{sample}/diamond_kegg/{sample}_diamond_kegg.tsv", # Or should it be named uniref-100?
    params: 
        query_cover = 85,
        subject_cover = 85, # 
        percent_id = 50, # 30 is probably fine for checkm2, but I feel like I'd rather have censored data than spurious results.
        evalue = "1e-05",
        blocksize = 2, # A value of 2 corresponds to running checkm2 in non-lowmem mode.
        database_path = DATABASES + "/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
    conda: "../envs/diamond.yaml" 
    resources:
        mem_mb = 20000, # Seems to use around 18G at max.
        runtime = "1h",
    benchmark: "{results_directory}/benchmarks/benchmark.diamond_kegg.{sample}.tsv"
    threads: 8
    shell: """

        # Inspired from https://github.com/chklovski/CheckM2/blob/319dae65f1c7f2fc1c0bb160d90ac3ba64ed9457/checkm2/diamond.py#L79
    
        # blastp: Align amino acid query sequences against a protein reference database

        diamond blastp \
            --outfmt 6 \
            --max-target-seqs 1 \
            --query {input.aminoacid:q}  \
            --out {output.tsv:q}  \
            --threads {threads}  \
            --db {params.database_path:q} \
            --query-cover {params.query_cover}  \
            --subject-cover {params.subject_cover}  \
            --id {params.percent_id}  \
            --evalue {params.evalue} \
            --block-size {params.blocksize}

    """
