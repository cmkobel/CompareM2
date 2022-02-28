Sys.setenv(LANG = "en") # Enforce english error messages on all systems.
message("loading tidyverse ..")
library(tidyverse)
message("load clusterProfiler..")
library(clusterProfiler) # Possible bug: I'm a bit worried that the download will fail on HPC systems. Not all HPC environments (i.e. singularity) allow internet access.. But: If it fails downloading, it can maybe use internal data?



args = commandArgs(trailingOnly = TRUE)
ko_file = args[1]
annotation_file = args[2]

# Development:
if (F) {
    ko_file = "~/assemblycomparator2/assets/ko"
    annotation_file = "~/assemblycomparator2/tests/E._faecium_plasmids/output_asscom2/collected_results/prokka_labelled.tsv"
}


message("these are the args: ", args)
message("this is ko_file:", ko_file)
message("this is annotation_file: ", annotation_file)




write("Reading ko asset", file = stderr())
kegg2gene = read_tsv(ko_file, col_names = c("ko_raw", "text")) %>% 
    separate(text, into = c("genes", "ko_text"), sep = ";") %>% # ignore warnings, it is likely because of ampersand encoding of weird characters we don't care about.
    mutate(gene = str_split(genes, ", ")) %>% 
    unnest(gene) %>% 
    mutate(ko = str_sub(ko_raw, 4)) %>% 
    
    # Finally, clean up:
    select(ko, gene, ko_text) %>% 
    identity()
    
    
write("Reading prokka labelled tsv", stderr())
annotations = read_tsv(annotation_file) %>% 
    rename(sample = last_col(), gene_raw = gene) %>%
    select(sample, everything()) %>% 
    filter(locus_tag != "locus_tag") %>% # Remove artifact
    mutate(length_bp = as.numeric(length_bp)) %>% 
    
    # Because Prokka marks different versions of the same gene, we need to remove appended underscores
    mutate(gene = str_remove(gene_raw, "_\\d+")) %>% 
    
    left_join(kegg2gene) %>% 
    
    drop_na(gene) %>% 
    identity()
    
    








# # Development: Explore the data a bit - Get an overview
# annotations %>% 
#     #drop_na(gene) %>% 
#     
#     ggplot(aes(gene, sample, fill = length_bp)) +
#     geom_tile() + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




# I. Perform GSEA on each sample
write("Initializing GSEA", stderr())
# Prepare the annotation data into a named list.
signatures = annotations %>% 
    #filter(sample == "Dallas_55") %>% # Develop; just one sample
    drop_na(ko) %>% 
    group_by(sample) %>% 
    summarize(genes = list(ko)) %>% 
    deframe()



# Development: Test that the enrichment works
# enrichKEGG(
#     gene = signatures$`Dallas_55`,
#     organism = "ko"
# ) %>% as_tibble()



# Loop over the samples and calculate the analysis independently for each of them.
sample_pathway_enrichment_analysis = tibble()
for (sample in names(signatures)) {
    write(paste("Performing enrichKEGG on sample:", sample), file = stderr())
    
    temp_anal = enrichKEGG(
        gene = signatures[[sample]],
        organism = "ko"
    ) %>% 
        as_tibble() %>% 
        mutate(sample = sample) %>% 
        rename_all(tolower) %>% 
        select(sample, everything())
    
    # Append the enrichment analysis
    sample_pathway_enrichment_analysis = sample_pathway_enrichment_analysis %>% 
        bind_rows(temp_anal)
    
}

# # Development: Let's get a crud viz
# # Make this figure in the report or wherever
# sample_pathway_enrichment_analysis %>% 
#     ggplot(aes(sample, description, fill = `p.adjust`)) +
#     geom_tile() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#     


# Write analysis to disk
write("Writing enrichment analysis to STDOUT", stderr())
sample_pathway_enrichment_analysis %>% 
    format_tsv() %>%
    cat()



# II. Perform set operations so as to showcase unique genes and their pathways in each sample
# I don't know what really makes sense. I think what I want to know is, which features are on the same plasmid. So I need to read the full gff file from prokka so I can put in the plasmid names. Maybe not something that should be visualized in the report, but it could be nice as a supplemental table.




