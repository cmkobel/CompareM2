
rule prokka:
    input: 
        metadata = "{results_directory}/metadata.tsv",
        assembly = "{results_directory}/samples/{sample}/{sample}.fna"
    output:
        gff = "{results_directory}/samples/{sample}/prokka/{sample}.gff",
        faa = "{results_directory}/samples/{sample}/prokka/{sample}.faa", # Used in dbcan, interproscan, diamond_kegg, motupan
        log = "{results_directory}/samples/{sample}/prokka/{sample}.log",
        tsv = "{results_directory}/samples/{sample}/prokka/{sample}.tsv",
        gff_nofasta = "{results_directory}/samples/{sample}/prokka/{sample}.gff_nofasta", # Not sure where I'm using this.
    conda: "../envs/prokka.yaml"
    benchmark: "{results_directory}/benchmarks/benchmark.prokka_individual.{sample}.tsv"
    resources:
        mem_mb = 8192,
    threads: 4
    shell: """
        
        prokka \
            --cpus {threads} \
            --force \
            --rfam \
            --compliant \
            --outdir {wildcards.results_directory}/samples/{wildcards.sample}/prokka \
            --prefix {wildcards.sample} {input.assembly:q} \
        | tee {output.log:q} 

        # I don't remember what I'm actually using this output for?
        # Remove fasta from gff and add sample label
        gff_fasta_start=$(grep --line-number --extended-regexp "^##FASTA" {output.gff:q} | cut -f1 -d:)
        head --lines $(($gff_fasta_start-1)) {output.gff:q} \
        > {output.gff_nofasta:q}

        {void_report}

    """
    