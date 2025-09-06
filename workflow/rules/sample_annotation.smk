
def get_annotation_results(wildcards):
    """ 
        Reads the config annotator parameter and returns the necessary 
        annotation results for that annotator. This could have also been an 
        anonymous rule but I think this is more snakemakistic. 
        This could also have been an action in the end of the annotaton rules (prokka and bakta), but then it would become complicated when the user wanted to run both, to possibly compare the output.
    """
    
    origin = df[df["sample"] == wildcards.sample]["origin"].tolist()[0]


    
    if origin == "local":
        
        raw = config["annotator"]
        parsed = str(raw).lower().strip()
        
        extra_message = "" # Extra message that can be added for debugging etc.
        
        # Instead of having this being set by the config["annotator"], it could be set by the corresponding line in the metadata table. Thus each sample could be either prokka, bakta or get_ncbi. It is fine to have the config set this, but then it should be copied into the metadata sheet first, and then taken from there.
        # The thing is that the annotator should not download, as we need the genome already when rule copy runs. So probably 
        if parsed == "bakta":
            annotator = "bakta"
        elif parsed == "prokka":
            annotator = "prokka"
        else: # base case, fall back.
            extra_message += " (default)"
            annotator = "bakta"
            
    elif origin == "ncbi":
        annotator = "ncbi_annotation"
            
        #print(f"Pipeline: Using the {annotator} annotator for sample \"{wildcards.sample}\"{extra_message}.")
        
        # Return should only contain one item, as I can't name them and I need to refer to a single one in rule annotate where I'm accessing its dirname().
    #print("Selected annotator for", wildcards.sample, ":", annotator)
    return [
            f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.gff",
            f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.faa",
            f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.log",
            f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.ffn",
            f"{wildcards.output_directory}/samples/{wildcards.sample}/{annotator}/{wildcards.sample}.tsv",
        ]
        

        
# rv = [
#     f"{wildcards.output_directory}/samples/{wildcards.sample}/ncbi/.renamed/{wildcards.sample}.gff",
#     f"{wildcards.output_directory}/samples/{wildcards.sample}/ncbi/.renamed/{wildcards.sample}.faa",
#     f"{wildcards.output_directory}/samples/{wildcards.sample}/ncbi/.renamed/{wildcards.sample}.log",
#     f"{wildcards.output_directory}/samples/{wildcards.sample}/ncbi/.renamed/{wildcards.sample}.ffn", # better check this one 
#     f"{wildcards.output_directory}/samples/{wildcards.sample}/ncbi/.renamed/{wildcards.sample}.tsv", # Not sure about this one. Maybe remove?
# ]
    

# There is a non-critical bug. When this rule is run after changing the annotator, the inputs from the old annotator are deleted? How does that work? Does snakemake keep track of changes in the job dag and remove old inputs?
# When the annotator choice is changed, the old files are deleted. they shouldn't though, so what I'm thinking is that these files should be linked as individual files instead of as a whole directory. Then the individual links can be deleted and new ones can be laid out. Problem solved, let's get cracking. Update: The problem with this solution is that it doesn't rerun the annotation step when the annotator is changed, but I guess that makes fine sense. Another update: Nu har jeg prøvet at skifte fra prokka til bakta, det ser ud til at den ikke laver nye links? Sandsynligvis fordi annotate slet ikke bliver kørt når man skriver --until bakta. Ja hvis man bare skriver --until annotate. After testing, I can confirm that it works well, but you will have to write --forcerun bakta if you want it to update. Conclusion: Now I've tried both approaches - both having a linked dir and individually linked files. I think a linked dir is cleaner, but the code is not because you need a lot of dirname commands, and you need to remove the (potential old) links in the end of bakta and prokka. The solution I have now, with individually linked files is simpler, so I think I'll stick with it.
rule annotate:
    input: 
        targets = get_annotation_results,
        #metadata = "{output_directory}/metadata.tsv", # Update linking when new genomes are input.
    output: # These are mostly the outputs that are used downstream.
        dir = directory("{output_directory}/samples/{sample}/.annotation/"),
        gff = "{output_directory}/samples/{sample}/.annotation/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/.annotation/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg, eggnog, amrfinderplus
        log = "{output_directory}/samples/{sample}/.annotation/{sample}.log",
        ffn = "{output_directory}/samples/{sample}/.annotation/{sample}.ffn",
        tsv = "{output_directory}/samples/{sample}/.annotation/{sample}.tsv",
    shell: """

        # Using a softlink means that the choice of annotator can be changed without loss of information. T
        # Though using this solution means that the user can mix up results, so that if the user switches between annotators while running the downstream tools, a mix of annotators can be used. This is very powerful in the sense that the user can use different annotators for different downstream tools, but it is also dangerous in the sense that it can create confusion around which tool is being used. Maybe there should be a log somewhere, or I should go back to having the directory being linked so that the "old" annotation results are removed everytime the annotator is switched.
        
        ln -sr {input} {output.dir}
        
    """
    
rule get_ncbi_annotation:
    input:
        assembly = "{output_directory}/samples/{sample}/{sample}.fna",
    output:
        gff_nofasta = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.gff_nofasta",
        gff = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.gff",
        
        faa = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.faa",
        log = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.log",
        ffn = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.ffn",
        tsv = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.tsv",
        gbk = "{output_directory}/samples/{sample}/ncbi_annotation/{sample}.gbk", 
    shell: """
    
        # Just a matter of linking the files to a meaningful directory. No annotation is performed, just that the files are put into a directory that looks a bit like the prokka and bakta directories.
        # Has to be copied rather than linked. So we don't link the same file for the non-nofasta gff.
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/genomic.gff {output.gff_nofasta}
        
        # Make a prokka-like gff containing the contigs after a ##FASTA line. 
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/genomic.gff {output.gff}
        echo "##FASTA" >> {output.gff}
        cat {input.assembly} >> {output.gff}
        
        
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/protein.faa {output.faa}
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/sequence_report.jsonl {output.log}
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/cds_from_genomic.fna {output.ffn}
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/genomic.gtf {output.tsv}
        cp {output_directory}/samples/{wildcards.sample}/ncbi/ncbi_dataset/data/genomic.gbff {output.gbk} # Not sure if it is problematic to have a "flat file". Let's see what happens.
        
    
        # TODO: Investigate what happens if a genbank or gtf file does not exist for a given accession. I simply don't know how comprehensive the ncbi db is.        


        
    
    """


rule prokka:
    input: 
        assembly = "{output_directory}/samples/{sample}/{sample}.fna",
    output:
        gff = "{output_directory}/samples/{sample}/prokka/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/prokka/{sample}.faa",
        ffn = "{output_directory}/samples/{sample}/prokka/{sample}.ffn",
        log = "{output_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{output_directory}/samples/{sample}/prokka/{sample}.tsv",
        gbk = "{output_directory}/samples/{sample}/prokka/{sample}.gbk",
        gff_nofasta = "{output_directory}/samples/{sample}/prokka/{sample}.gff_nofasta", # Might come in handy.
    params:
        passthrough_parameters = passthrough_parameter_unpack("prokka")
    conda: "../envs/prokka.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.prokka_sample.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
    
        # Collect version number.
        prokka --version > "$(dirname {output.gff})/.software_version.txt"
        
        prokka \
            --cpus {threads} \
            --force \
            {params.passthrough_parameters} \
            --outdir {wildcards.output_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input.assembly:q} \
        | tee {output.log:q} 
        
        

        # Remove fasta from gff and add sample label
        sed '/^##FASTA/Q' {output.gff:q} > {output.gff_nofasta:q}

        {void_report}

    """
    

rule bakta:
    input: 
        database_representative = DATABASES + "/bakta/comparem2_bakta_database_representative.flag",
        assembly = "{output_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{output_directory}/samples/{sample}/bakta/{sample}.gff",
        faa = "{output_directory}/samples/{sample}/bakta/{sample}.faa",
        tsv = "{output_directory}/samples/{sample}/bakta/{sample}.tsv",
        log = "{output_directory}/samples/{sample}/bakta/{sample}.log",
        ffn = "{output_directory}/samples/{sample}/bakta/{sample}.ffn",
        #gff_generic = "{output_directory}/samples/{sample}/annotation/{sample}.gff3",
        # no gbk file?
    params:
        DATABASES = DATABASES,
        passthrough_parameters = passthrough_parameter_unpack("bakta")
    conda: "../envs/bakta.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.bakta_sample.{sample}.tsv"
    resources:
        mem_mb = 16384,
    threads: 4
    shell: """
                        
        # Collect version number.
        bakta --version > "$(dirname {output.gff})/.software_version.txt"
        
        # Collect database version.
        echo -e "$(date -Iseconds)\t$(dirname {input.database_representative})" > "$(dirname {output.gff})/.database_version.txt"
        
        bakta \
            --db {params.DATABASES}/bakta/db \
            --output $(dirname {output.gff}) \
            --threads {threads} \
            --force \
            {params.passthrough_parameters} \
            {input.assembly}
            
        # Optimize compatibility between prokka and bakta to make fluent use of either result easy.
        cp {output.gff}3 {output.gff}
        
        {void_report}

    """


# rule get_ncbi:
#     output:
#         gff = "{output_directory}/samples/{sample}/bakta/{sample}.gff",
#         faa = "{output_directory}/samples/{sample}/bakta/{sample}.faa",
#         tsv = "{output_directory}/samples/{sample}/bakta/{sample}.tsv",
#         log = "{output_directory}/samples/{sample}/bakta/{sample}.log",
#         ffn = "{output_directory}/samples/{sample}/bakta/{sample}.ffn",
#     conda: "../envs/ncbi-datasets.yaml"
#     benchmark: "{output_directory}/benchmarks/benchmark.get_ncbi_sample.{sample}.tsv"
#     shell: """
    
#         datasets \
#             download genome \
#             accession GCF_900186865.1 \
#             --include gff3,rna,cds,protein,genome,seq-report
    
#     """
