
# By default it runs TIGRFAM, Hamap, Pfam
rule interproscan:
    input: 
        aminoacid = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
        #metadata = "{output_directory}/metadata.tsv",
        # No external database is needed.
    output:
        tsv = "{output_directory}/samples/{sample}/interproscan/{sample}_interproscan.tsv",
    params:
        file_base = "{output_directory}/samples/{sample}/interproscan/{sample}_interproscan",
        passthrough_parameters = passthrough_parameter_unpack("interproscan")
    conda: "../envs/interproscan.yaml" # Not sure if it should be called by a version number?
    benchmark: "{output_directory}/benchmarks/benchmark.interproscan.{sample}.tsv"
    threads: 16
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
        aminoacid = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
        database_representative = DATABASES + "/dbcan/comparem2_dbcan_database_representative.flag"
        #metadata = "{output_directory}/metadata.tsv",
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


# Why am I using both results and output directory here?

# aka eggnog-mapper
rule eggnog:
    input: 
        database_representative = DATABASES + "/eggnog/comparem2_eggnog_database_representative.flag",
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg, eggnog
        #assembly = "{output_directory}/samples/{sample}/{sample}.fna",
        #metadata = "{output_directory}/metadata.tsv",
    output:
        gff = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.decorated.gff", # Looks like not genepred, but decorated is produced. Is that because of the settings or version changes?
        hits = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.hits",
        orthologs = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.seed_orthologs",
        tsv = "{output_directory}/samples/{sample}/eggnog/{sample}.emapper.annotations",
        # 2.1.12: {sample}.emapper.decorated.gff, {sample}.emapper.hits, {sample}.emapper.seed_orthologs, {sample}.emapper.annotations,
    params:
        passthrough_parameters = passthrough_parameter_unpack("eggnog")
    conda: "../envs/eggnog.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.eggnog_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 16 # At least for diamond, parellel computation is efficient.
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
            -i {input.faa:q} \
            --cpu {threads} \
            --itype proteins \
            --override \
            --temp_dir $TMPDIR \
            --decorate_gff yes \
            {params.passthrough_parameters}
            
        # Originally I wanted eggnog to be an alternative to prokka or bakta, but since it doesn't produce a .faa file there is no point in using like so. So now there is no point in re-calling genes with prodigal at genome state - I might as well reuse the protein .faa from the annotation tool instead. This should be faster
            
        # touch {output} # Just to check what comes out.

        {void_report}

    """
    

# Finally, gapseq is on bioconda, and seems to work. But I now don't understand why I did the rule order thing - Is that really necessary? Seams unnecessarily complicated.
rule gapseq_find:
    input: 
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa"
        #metadata = "{output_directory}/metadata.tsv",
    output:
        dir = directory("{output_directory}/samples/{sample}/gapseq"),
        pathways = "{output_directory}/samples/{sample}/gapseq/{sample}-all-Pathways.tbl",
        reactions = "{output_directory}/samples/{sample}/gapseq/{sample}-all-Reactions.tbl",
        transporter = "{output_directory}/samples/{sample}/gapseq/{sample}-Transporter.tbl",
    conda: "../envs/gapseq.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.gapseq_find.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 2
    shell: """
    
        # Collect version number.
        gapseq -v | awk '{{ print "gapseq", $0 }}' > "$(dirname {output.pathways})/.software_version.txt"
        
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

        {void_report}

    """
    
    
rule gapseq_fill: # Continuation on gapseq_find results.
    input: 
        assembly = "{output_directory}/samples/{sample}/{sample}.fna",
        pathways = "{output_directory}/samples/{sample}/gapseq/{sample}-all-Pathways.tbl",
        reactions = "{output_directory}/samples/{sample}/gapseq/{sample}-all-Reactions.tbl",
        transporter = "{output_directory}/samples/{sample}/gapseq/{sample}-Transporter.tbl",
        #metadata = "{output_directory}/metadata.tsv",
    output:
        rxnWeights = "{output_directory}/samples/{sample}/gapseq/{sample}-rxnWeights.RDS", 
        rxnXgenes = "{output_directory}/samples/{sample}/gapseq/{sample}-rxnXgenes.RDS",
        draft = "{output_directory}/samples/{sample}/gapseq/{sample}-draft.RDS",
        draft_xml = "{output_directory}/samples/{sample}/gapseq/{sample}-draft.xml",
    
        medium = "{output_directory}/samples/{sample}/gapseq/{sample}-medium.csv",
    
        #filled = "{output_directory}/samples/{sample}/gapseq/{sample}.tbl",
        filled_xml = "{output_directory}/samples/{sample}/gapseq/{sample}.xml"
    params:
        dir = "{output_directory}/samples/{sample}/gapseq",
        passthrough_parameters_draft = passthrough_parameter_unpack("gapseq_fill_draft"),
        passthrough_parameters_medium = passthrough_parameter_unpack("gapseq_fill_medium"),
        passthrough_parameters_fill = passthrough_parameter_unpack("gapseq_fill_fill"),
        gapseq_medium = gapseq_medium
    conda: "../envs/gapseq.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.gapseq.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 1 # Gapseq is very badly optimized, so mostly runs with 4% CPU load on a single core anyway. 
    shell: """
    
        echo "1) Drafting ..."

        # Produces *-rxnWeights.RDS, *-rxnXgenes.RDS, *-draft.RDS, *-draft.xml
        gapseq draft \
            -r {input.reactions:q} \
            -t {input.transporter:q} \
            -p {input.pathways:q} \
            -c {input.assembly:q} \
            -f {params.dir:q} \
            {params.passthrough_parameters_draft}
            

        echo "2) Preparing medium definition ..."
                
        if [ -f "{params.gapseq_medium}" ]; then
            echo "Medium: Predefined path exists. Copying '{params.gapseq_medium}'"
            cat "{params.gapseq_medium}" > {output.medium}
            echo -e "gapseq\tmedium {params.gapseq_medium}" >> "$(dirname {output.medium})/.software_version.txt"
        else
            echo "Medium: Predicting growth medium ..."
            gapseq medium -m {output.draft} -p {input.pathways} -o {output.medium} {params.passthrough_parameters_medium}
            echo -e "gapseq:\tmedium predicted de novo using gapseq medium" > "$(dirname {output.medium})/.software_version.txt"

        fi
        

        echo "3) Filling ..."    
        
        # Produces *.RDS and *.xml
        gapseq fill \
            -m {output.draft:q} \
            -c {output.rxnWeights:q} \
            -g {output.rxnXgenes:q} \
            -f {params.dir:q} \
            -n {output.medium} \
            {params.passthrough_parameters_fill}
            

        # {void_report} How to present this?

    """
    
# This is a simple idea. Use the dram standard report as basis. Enumerate the same pathways/functions and show, with metabolic rates, which are present based on a specific medium. Could be relevant for streptococcus.
# rule gapseq_survey:
#     input: filled_xml = "{output_directory}/samples/{sample}/gapseq/{sample}.xml"
#     output: survey_tsv = "{output_directory}/samples/{sample}/gapseq/{sample}_survey.tsv
#     conda: "../envs/r-cobrar.yaml"
#     benchmark: "{output_directory}/benchmarks/benchmark.gapseq.{sample}.tsv"
#     shell: """ # partly copied from batch_advanced_annotation::rule kegg_pathway
    
#         # Collect version number.
#         R -s -q -e "library(cobrar); sessionInfo()"  | grep -P "R version|clusterProfiler" > "$(dirname {output.pathway_enrichment})/.software_version.txt"

#         Rscript {params.script:q} \
#             {input.kegg_asset:q} \
#             {params.output_dir:q} \
#             {input.eggnog_annotations:q}
            
#         {void_report}

#     """
    
    

    
    
    

rule antismash:
    input: 
        database_representative = DATABASES + "/antismash/comparem2_antismash_database_representative.flag",
        gbk = "{output_directory}/samples/{sample}/prokka/{sample}.gbk", #culprit for always running prokka? Does bakta produce a gbk?
        #metadata = "{output_directory}/metadata.tsv",
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
        rm -rf $(dirname {output.json})/* || echo "Continuing ..."
        
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

        # {void_report} Not yet implemented.

    """
    
    
rule carveme: # Gapseq is too slow.
    input: 
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa",
    output:
        xml = "{output_directory}/samples/{sample}/carveme/{sample}_carveme.xml",
    conda: "../envs/carveme.yaml"
    params:
        passthrough_parameters = passthrough_parameter_unpack("carveme"),
    benchmark: "{output_directory}/benchmarks/benchmark.carveme_sample.{sample}.tsv"
    shell: """
    
        carve --help || echo OK
        
        # Build a model and gapfil it proritizing the reactions selected for gap-filling based on genetic evidence.
        carve {input.faa} \
            -o {output.xml} \
            -g M9 \
            -i M9 \
            {params.passthrough_parameters}
        
        
        # Microbial communities (another job?)
    
    
    
    """