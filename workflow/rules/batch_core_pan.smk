
    
        
        
def get_mem_panaroo(wildcards, attempt): 
    return [32000, 64000, 128000][attempt-1]



#print(passthrough_parameter_unpack("none")) # DEBUG

rule panaroo: # Checkpoint, because some rules i.e. fasttree, iqtree, snp-dists should only run if the core genome is non-empty.
    input: 
        metadata = "{output_directory}/metadata.tsv",
        gff = expand("{output_directory}/samples/{sample}/.annotation/{sample}.gff", sample = df["sample"], output_directory = output_directory),
    output:
        summary = "{output_directory}/panaroo/summary_statistics.txt", # Todo: don't rely on this file as it is not produced when the core is empty. Instead, the core should be inferred from the gene_absence_presence file.
        presence = "{output_directory}/panaroo/gene_presence_absence.Rtab",
        #alignment = "{output_directory}/panaroo/core_gene_alignment.aln",
        #analyses = ["{output_directory}/panaroo/summary_statistics.txt", "{output_directory}/panaroo/core_gene_alignment.aln", "{output_directory}/panaroo/gene_presence_absence.csv"]
    params:        
        script_summary = base_variable + "/workflow/scripts/panaroo_generate_summary_stats.py",
        passthrough_parameters = passthrough_parameter_unpack("panaroo")
    benchmark: "{output_directory}/benchmarks/benchmark.panaroo.tsv"
    threads: 16
    #retries: 2
    resources:
        #mem_mb = 32768,
        mem_mb = get_mem_panaroo,
        runtime = "24h",
    conda: "../envs/panaroo.yaml" 
    shell: """
    
        # Collect version number.
        panaroo --version > "$(dirname {output.summary})/.software_version.txt"
            
        panaroo \
            -o {wildcards.output_directory}/panaroo \
            -t {threads} \
            {params.passthrough_parameters} \
            -i {input.gff:q}

        # The summary stat table is very handy to quickly check the size of the core genome. Unfortunately it is not generated when the core genome is absent, so here comparem2 is reproducing it from the presence file which is generated in all cases.
        # A better solution might actually be to change panaroo so it produces this file nonetheless. 
        #{params.script_summary} {output.presence:q} > {output.summary:q}
        # Nevermind, turns out that it is the alignment file that isn't created, which makes more sense. # Now I just need to find a way of returning the alignment file when that exists.
        
        {void_report}

    """
    

# Instead of passing a missing file into the dependents of the core genome, this rule fails and stops the downstream jobs. Note that this is not a "real" snakemake-checkpoint, as those cannot support conditional running whatsoever.
# A possible todo is to make this into a python script that parses the summary and actually prints the size of the core genome. But alas, another day.
rule core_genome_checkpoint:
    input: 
        summary = "{output_directory}/panaroo/summary_statistics.txt"
    output: 
        aln = "{output_directory}/panaroo/core_gene_alignment_verified.aln"
    shell: """
    
        if [ -f "{output_directory}/panaroo/core_gene_alignment.aln" ]; then
            echo "Pipeline: Core genome exists."
            cp "{output_directory}/panaroo/core_gene_alignment.aln" {output.aln}
            ls -lh {output.aln} # show size
        else 
            echo "Pipeline: Warning: Core genome does not exist: Downstream analyses will be cancelled."
            exit 1
        fi    
        
    """


rule snp_dists:
    input: 
        metadata = "{output_directory}/metadata.tsv",
        aln = "{output_directory}/panaroo/core_gene_alignment_verified.aln",
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
