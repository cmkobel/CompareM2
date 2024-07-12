
    
        
        
def get_mem_panaroo(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]



#print(passthrough_parameter_unpack("none")) # DEBUG

# Maybe it shouldn't be panaroo which is the checkpoint, but rather every rule that uses the core genome, like rule snp_dists? Let me try!
checkpoint panaroo: # Checkpoint, because some rules i.e. fasttree, iqtree, snp-dists should only run if the core genome is non-empty.
    input: 
        metadata = "{output_directory}/metadata.tsv",
        gff = expand("{output_directory}/samples/{sample}/prokka/{sample}.gff", sample = df["sample"], output_directory = output_directory),
    output:
        summary = "{output_directory}/panaroo/summary_statistics.txt",
        presence = "{output_directory}/panaroo/gene_presence_absence.csv",
        alignment = "{output_directory}/panaroo/core_gene_alignment.aln",
        #analyses = ["{output_directory}/panaroo/summary_statistics.txt", "{output_directory}/panaroo/core_gene_alignment.aln", "{output_directory}/panaroo/gene_presence_absence.csv"]
    params:        
        passthrough_parameters = passthrough_parameter_unpack("panaroo")
    benchmark: "{output_directory}/benchmarks/benchmark.panaroo.tsv"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_panaroo,
        runtime = "24h",
    conda: "../envs/panaroo.yaml" # deactivated 
    shell: """
    
        # Collect version number.
        panaroo --version > "$(dirname {output.summary})/.software_version.txt"
            
        panaroo \
            -o {wildcards.output_directory}/panaroo \
            -t {threads} \
            {params.passthrough_parameters} \
            -i {input.gff:q}
        
        {void_report}

    """
    

def core_genome_if_exists(wildcards): 
    
    summary_file = checkpoints.panaroo.get(**wildcards).output["summary"]
    alignment_file = checkpoints.panaroo.get(**wildcards).output["alignment"]
    
    try:   
        with summary_file.open() as f:
        #with open(file) as f: 

            for line in f:
                line = line.strip()
                #print(f"line: {line}")
                
                # Find the line that states the size of the core genome.
                if "Core genes" in line:
                    #print(f"Pipeline: Parsing Panaroo core genome from line \"{line}\"")
                    splitted = line.split("\t")
                    #print(splitted)
                    core_genome_size = int(splitted[-1])
                    
                    #print(f"Pipeline: Parsed core genome size is: {repr(core_genome_size)}")
                    break
                
            # Return alignment file when there is a core genome.
            if core_genome_size > 0:
                print(f"Pipeline: Core genome exists ({core_genome_size} genes).")
                return [alignment_file]
                
            # Otherwise, return an empty list.
            else:
                print(f"Pipeline: Core genome is empty ({core_genome_size} genes).")
                return list()
                
                
    # In case the file cannot be parsed, we're assuming that there is no core genome. Hence return empty list.
    except Exception as e: 
        print(f"Pipeline Error: {e}") # Show errors in parsing the panaroo summary, but don't fail.
        print("Pipeline: Core genome size not known?")
        return list() # Empty list
        



rule compute_snp_dists:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        aln = core_genome_if_exists,
    output: "{output_directory}/snp-dists/snp-dists.tsv"
    conda: "../envs/snp-dists.yaml"
    benchmark: "{output_directory}/benchmarks/benchmark.snp_dists.tsv"
    threads: 4
    shell: """
    
        # Collect version number.
        snp-dists -v > "$(dirname {output})/.software_version.txt"

        snp-dists \
            -j {threads} \
            {input.aln:q} > {output:q}

        {void_report}
        
    """

# This is a dummy rule that forces recomputation of the DAG after the core genome alignment has been produced. It solves a possibly known error in snakemake which I could possibly mitigate by updating snakemake (TODO). But since I'm already testing a lot of other features, I cannot take the burden of also changing snakemake version so this is the workaround for now.
rule snp_dists:
    input: "{output_directory}/snp-dists/snp-dists.tsv"
    output: touch("{output_directory}/snp-dists/.done.flag")
    shell: """
        exit 0
    """
