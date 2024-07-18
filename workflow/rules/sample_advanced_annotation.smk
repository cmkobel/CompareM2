
# By default it runs TIGRFAM, Hamap, Pfam
rule interproscan:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        aminoacid = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
        # No external database is needed.
    output:
        tsv = "{output_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv",
    params:
        file_base = "{output_directory}/samples/{sample}/interproscan/{sample}_interproscan",
        passthrough_parameters = passthrough_parameter_unpack("interproscan")
    conda: "../envs/interproscan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{output_directory}/benchmarks/benchmark.interproscan.{sample}.tsv"
    threads: 8
    resources: 
        mem_mb = 8000
    shell: """
    
        # Collect version number.
        interproscan.sh --version | grep version > "$(dirname {output.tsv})/.software_version.txt"

        # https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html#command-line-options
        interproscan.sh \
            --cpu {threads} \
            --output-file-base {params.file_base:q} \
            --disable-precalc \
            --formats TSV \
            --iprlookup \
            {params.passthrough_parameters} \
            --seqtype p \
            --tempdir {resources.tmpdir} \
            --input {input.aminoacid:q}

    """


rule dbcan: # I can't decide whether this rule should really be called "run_dbcan", since that is the name of the software.
    input: 
        metadata = "{output_directory}/metadata.tsv",
        aminoacid = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
        database_representative = DATABASES + "/dbcan/comparem2_dbcan_database_representative.flag"
    output: 
        overview_table = "{output_directory}/samples/{sample}/dbcan/overview.txt",
        diamond_table = "{output_directory}/samples/{sample}/dbcan/diamond.out"
    params: 
        out_dir = "{output_directory}/samples/{sample}/dbcan"
    conda: "../envs/dbcan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{output_directory}/benchmarks/benchmark.dbcan.{sample}.tsv"
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




# aka eggnog-mapper
rule eggnog:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        database_representative = DATABASES + "/eggnog/comparem2_eggnog_database_representative.flag",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        hits = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.hits",
        orthologs = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.seed_orthologs",
        tsv = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.annotations",
        ffn = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.genepred.fasta",
        
        gff = "{results_directory}/samples/{sample}/eggnog/{sample}.emapper.gff", # Why is it sometimes called emapper.genepred.gff?
    params:
        passthrough_parameters = passthrough_parameter_unpack("eggnog")
    conda: "../envs/eggnog.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.eggnog_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 8 # Not sure if the underlying tools are capable of doing lots of parallel computation.
    shell: """
        
        # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#basic-usage
        
        # Collect version number.
        emapper.py \
            --data_dir $(dirname {input.database_representative}) \
            --version > "$(dirname {output.gff})/.software_version.txt"
            
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.gff})/.database_version.txt"
        
        emapper.py \
            --data_dir $(dirname {input.database_representative}) \
            --output_dir "$(dirname {output.gff})/" \
            -o "{wildcards.sample}" \
            -i {input.assembly:q} \
            --cpu {threads} \
            --itype genome \
            --override \
            --temp_dir $TMPDIR \
            --decorate_gff yes \
            {params.passthrough_parameters}
            
        touch {output} # Just to check what comes out.

        {void_report}

    """
    


rule gapseq_find:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa"
    output:
        dir = directory("{output_directory}/samples/{sample}/gapseq"),
        pathways = "{output_directory}/samples/{sample}/gapseq/{sample}-Pathways.tbl",
        reactions = "{output_directory}/samples/{sample}/gapseq/{sample}-Reactions.tbl",
        transporter = "{output_directory}/samples/{sample}/gapseq/{sample}-Transporter.tbl",
        flag = "{output_directory}/samples/{sample}/gapseq/gapseq_done.flag",
    conda: "../envs/gapseq.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.gapseq_find_sample.{sample}.tsv"
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
        metadata = "{output_directory}/metadata.tsv",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna",
        pathways = "{output_directory}/samples/{sample}/gapseq/{sample}_pathways.tsv",
    output:
        rxnWeights = "{output_directory}/samples/{sample}/gapseq/{sample}-rxnWeights.tbl",
        rxnXgenes = "{output_directory}/samples/{sample}/gapseq/{sample}-rxnXgenes.tbl",
        draft = "{output_directory}/samples/{sample}/gapseq/{sample}-draft.tbl",
        draft_xml = "{output_directory}/samples/{sample}/gapseq/{sample}-draft_xml.xml",

        filled = "{output_directory}/samples/{sample}/gapseq/{sample}.tbl",
        filled_xml = "{output_directory}/samples/{sample}/gapseq/{sample}.xml",
        
        flag = "{output_directory}/samples/{sample}/gapseq/gapseq_done.flag",
    params:
        dir = "{output_directory}/samples/{sample}/gapseq",

    conda: "../envs/gapseq.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.gapseq_find_sample.{sample}.tsv"
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
        metadata = "{output_directory}/metadata.tsv",
        database_representative = DATABASES + "/antismash/comparem2_antismash_database_representative.flag",
        gbk = "{output_directory}/samples/{sample}/prokka/{sample}.gbk",
    output:
        json = "{output_directory}/samples/{sample}/antismash/{sample}.json",
    params:
        DATABASES = DATABASES,
        dir = "{output_directory}/samples/{sample}/antismash",
    conda: "../envs/antismash.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.antismash_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 8
    shell: """
    
        # Clean directory. Antismash will fail if previous files exist.
        rm $(dirname {output.json:q})/* || echo no files
        
        # Collect version number.
        antismash --version > "$(dirname {output.json})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.json})/.database_version.txt"
    
        antismash \
            --taxon bacteria \
            --cpus {threads} \
            --output-dir {params.dir:q} \
            --output-basename {wildcards.sample:q} \
            --databases "{params.DATABASES}/antismash" \
            {input.gbk:q}

        {void_report}

    """
    