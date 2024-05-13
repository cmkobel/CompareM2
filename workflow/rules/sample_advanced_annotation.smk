
# By default it runs TIGRFAM, Hamap, Pfam
rule interproscan:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/.annotation/{sample}.faa",
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
    
        # Collect version number.
        interproscan.sh --version | grep version > "$(dirname {output.tsv})/.software_version.txt"

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



# Use the same database as checkm2, but run on amino-acid files instead of dna.
# This will be used for a subsequent pathway enrichment analysis.
# Idea for speed up: concatenate all genomes together first, like they do in checkm2. Then we only need to load the database once.
rule diamond_kegg: # or uniref_ko?
    input: 
        metadata = "{results_directory}/metadata.tsv", # For the report
        aminoacid = "{results_directory}/samples/{sample}/.annotation/{sample}.faa", 
        #database_representative = base_variable + "/databases/checkm2/ac2_checkm2_database_representative.flag",
        database_representative = DATABASES + "/checkm2/ac2_checkm2_database_representative.flag",
        
    output: # 
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
        
        # Collect version number.
        diamond version > "$(dirname {output.tsv})/.software_version.txt"

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

   

rule dbcan: # I can't decide whether this rule should really be called "run_dbcan", since that is the name of the software.
    input: 
        metadata = "{results_directory}/metadata.tsv",
        aminoacid = "{results_directory}/samples/{sample}/.annotation/{sample}.faa",
        database_representative = DATABASES + "/dbcan/ac2_dbcan_database_representative.flag"
    output: 
        overview_table = "{results_directory}/samples/{sample}/dbcan/overview.txt",
        diamond_table = "{results_directory}/samples/{sample}/dbcan/diamond.out"
    params: 
        out_dir = "{results_directory}/samples/{sample}/dbcan"
    conda: "../envs/dbcan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{results_directory}/benchmarks/benchmark.dbcan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """
    
        # Collect version number.
        conda list | grep dbcan > "$(dirname {output.overview_table})/.software_version.txt" # Todo, test on docker.
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.overview_table})/.database_version.txt"

        # It seems to be necessary to set all the cpu thread counts manually.
        export HMMER_NCPU={threads}

        run_dbcan \
            --dbcan_thread {threads} \
            --dia_cpu {threads} \
            --hmm_cpu {threads} \
            --tf_cpu {threads} \
            --stp_cpu {threads} \
            --db_dir $(dirname {input.database_representative:q}) \
            --out_dir {params.out_dir:q} \
            {input.aminoacid:q} \
            protein 

    """






rule gapseq_find:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        faa = "{results_directory}/samples/{sample}/.annotation/{sample}.faa"
    output:
        dir = directory("{results_directory}/samples/{sample}/gapseq"),
        pathways = "{results_directory}/samples/{sample}/gapseq/{sample}-Pathways.tbl",
        reactions = "{results_directory}/samples/{sample}/gapseq/{sample}-Reactions.tbl",
        transporter = "{results_directory}/samples/{sample}/gapseq/{sample}-Transporter.tbl",
        flag = "{results_directory}/samples/{sample}/gapseq/gapseq_done.flag",
    conda: "../envs/gapseq.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gapseq_find_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
    
    
        # Collect version number.
        # Todo
        
        # -K is only for multiple sequence alignments
        # -O is "offline mode"
        
        # Produces *-Pathways.tbl, *-Reactions.tbl
        gapseq find \
            -M prot \
            -f {output.dir} \
            -K {threads} \
            -O \
            -p all \
            {input.faa:q} 
            
        # Produces *-Transporter.tbl
        gapseq find-transport \
            -M prot \
            -f {output.dir} \
            -K {threads} \
            {input.faa:q}
            
            
        touch {output.flag} # To have a common output with rule gapseq for the ruleorder
            
        touch {output} # DEBUG
    

        {void_report}

    """
    
    
rule gapseq: # Continuation on gapseq_find results.
    input: 
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna",
        pathways = "{results_directory}/samples/{sample}/gapseq/{sample}_pathways.tsv",
    output:
        rxnWeights = "{results_directory}/samples/{sample}/gapseq/{sample}-rxnWeights.tbl",
        rxnXgenes = "{results_directory}/samples/{sample}/gapseq/{sample}-rxnXgenes.tbl",
        draft = "{results_directory}/samples/{sample}/gapseq/{sample}-draft.tbl",
        draft_xml = "{results_directory}/samples/{sample}/gapseq/{sample}-draft_xml.xml",

        filled = "{results_directory}/samples/{sample}/gapseq/{sample}.tbl",
        filled_xml = "{results_directory}/samples/{sample}/gapseq/{sample}.xml",
        
        flag = "{results_directory}/samples/{sample}/gapseq/gapseq_done.flag",
    params:
        dir = "{results_directory}/samples/{sample}/gapseq",

    conda: "../envs/gapseq.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.gapseq_find_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """

        # Produces *-rxnWeights.RDS, *-rxnXgenes.RDS, *-draft.RDS, *-draft.xml
        gapseq draft \
            -r toy/myb71-all-Reactions.tbl \
            -t toy/myb71-Transporter.tbl \
            -p toy/myb71-all-Pathways.tbl \
            -c {input:assembly:q}
        
        # Produces *.RDS and *.xml
        gapseq fill \
            -m toy/myb71-draft.RDS \
            -c toy/myb71-rxnWeights.RDS \
            -g toy/myb71-rxnXgenes.RDS \
            -n dat/media/TSBmed.csv
            
            
            
        touch {output.flag} # To have a common output with rule gapseq for the ruleorder
        
        touch {output} # DEBUG

        {void_report}

    """
    
    

rule antismash:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/antismash/ac2_antismash_database_representative.flag",
        fna = "{results_directory}/samples/{sample}/{sample}.fna",
    output:
        json = "{results_directory}/samples/{sample}/antismash/{sample}.json",
    params:
        DATABASES = DATABASES,
        dir = "{results_directory}/samples/{sample}/antismash",
    conda: "../envs/antismash.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.antismash_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 8
    shell: """
    
        # Collect version number.
        antismash --version > "$(dirname {output.json})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.json})/.database_version.txt"
    
        antismash \
            -t bacteria \
            -c {threads} \
            --output-dir {params.dir:q} \
            --output-basename {wildcards.sample:q} \
            --databases "{params.DATABASES}/antismash" \
            {input.fna:q}

        {void_report}

    """
    