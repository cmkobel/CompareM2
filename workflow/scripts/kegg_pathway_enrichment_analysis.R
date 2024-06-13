library(tidyverse)
library(clusterProfiler)
library(jsonlite)


# Present arguments for debugging
# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)


input_kegg_asset <- args[1]
output_path <- args[2]
input_diamond <- args[3:(length(args))] # Remaining arguments.



# For development
if (F) {
    # Work macbook
    input_kegg_asset <- "~/comparem2/resources/ko00001.json"
    output_path <- paste0(getwd(), "/outt")
    # input_diamond = "results_comparem2/samples/Treatment_LowE.metabat.242/diamond_kegg/Treatment_LowE.metabat.242_diamond_kegg.tsv"
    input_diamond <- "samples/116_2/eggnog/116_2.emapper.annotations"
}


# message(paste(args, collapse = "\n"))
message("Command line arguments: (", length(args), "):")
message("  input_kegg_asset: ", input_kegg_asset) # Path to json downloaded from kegg.jp
message("  output_path: ", output_path) # A path where to save analysis results.
message("  input_diamond (", length(input_diamond), "): ", paste(input_diamond, collapse = ", ")) # Table from diamond containing calls on the uniref100-KO database.
message("")



## Parse KEGG json

kegg_parse_json <- function(input_json) {
    kegg <- jsonlite::fromJSON(input_json, flatten = F)

    unnest <- kegg$children %>%
        unnest(children, names_repair = "universal") %>% # AB 57
        unnest(children, names_repair = "universal") %>% # C 551
        unnest(children, names_repair = "universal") %>% # D 59868

        rename(
            class = 1, # E.g. "09100 Metabolism"                # aka. class?
            group = 2, # E.g. "09101 Carbohydrate metabolism"   # aka. group?
            pathway = 3,
            ortholog = 4
        ) %>%
        separate_wider_regex(pathway, # E.g. "00010 Glycolysis \/ Gluconeogenesis [PATH:ko00010]"
            c(
                pathway = "^[^\\[]+", # Anything that isn't a [
                pathway_id = "\\[.+\\]$"
            ),
            too_few = "align_start"
        ) %>%
        separate_wider_regex(ortholog, # E.g. "K00844  HK; hexokinase [EC:2.7.1.1]"
            c( # "name":"K04022  eutG; alcohol dehydrogenase"
                ortholog_ko = "^K[0-9]+",
                "  ",
                ortholog_abbreviation = "[^;]+", # Anything that isn't a ;
                "; ",
                ortholog_name = "[^\\[]+", # Anything that isn't a [
                " \\[",
                ortholog_ec = ".+",
                "\\]$"
            ),
            too_few = "align_start"
        ) %>%
        mutate(ortholog_ec = str_trim(ortholog_ec)) %>%
        mutate(illegal_ortholog = is.na(ortholog_ko))

    # Check and report on the number of orthologs with illegal names.
    illegals <- unnest %>% filter(illegal_ortholog)
    if (nrow(illegals) > 0) {
        message(paste("\nWarning:", nrow(illegals), "ortholog(s) have an illegal name."))
        message(illegals$ortholog %>% head())
    }

    # Return everything that isn't illegal.
    unnest %>%
        filter(!illegal_ortholog) %>%
        select(-illegal_ortholog)
}



# Turns out you can just download this json: https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir=
kegg_data <- kegg_parse_json(input_kegg_asset)

# Debug
kegg_data %>% glimpse()

# Save the kegg data in case the user wants to investigate the pathway hierarchy manually.
kegg_data %>%
    write_tsv(paste0(output_path, "/kegg_data.tsv"))


# TODO:
# Make a sanity check that makes sure that there aren't too many NAs or something.


## Parse diamond results that link the called prodigal-called proteins to uniref100-KO
eggnog_raw <- tibble(File = Sys.glob(input_diamond)) %>% #
    extract(File, "sample", "/eggnog/(.+)\\.emapper\\.annotations", remove = F) %>%
    mutate(
        tabulation = lapply(
            File,
            read_tsv,
            col_names = c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"),
            comment = "#",
            na = "-"
        )
    ) %>%
    select(-File) %>%
    unnest(tabulation)





eggnog_raw %>%
    count(sample) %>%
    glimpse()


df_ko <- eggnog_raw %>%
    select(sample, KEGG_ko) %>%
    drop_na(KEGG_ko) %>%
    mutate(koid = str_split(KEGG_ko, ",")) %>%
    unnest(koid) %>%
    mutate(koid = str_remove(koid, "^ko:")) %>%
    select(sample, koid)






## Universal PEA
# https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html#
# Prepare genelist
# https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelis

# enricher


# term2gene is the geneset
term2gene <- kegg_data %>%
    select(term = pathway, gene = ortholog_ko)


# Then for each sample

analyses <- tibble() # For collecting the results.
for (sample in df_ko %>%
    group_by(sample) %>%
    group_split()) {
    sample_name <- sample %>%
        pull(sample) %>%
        unique()
    message("Running clusterProfiler::enricher on ", sample_name)

    # gene is the signature
    gene <- sample %>%
        pull(koid)

    # This is where we call the enrichment analysis.
    analysis <- clusterProfiler::enricher(
        gene,
        TERM2GENE = term2gene
    ) %>%
        as_tibble() %>%
        mutate(sample = sample_name) %>%
        relocate(sample) %>% # Put sample column in front.
        select(-Description) %>% # Redundant with ID
        rename(pathway = ID)


    # Collect into one big table.
    analyses <- bind_rows(
        analyses, # Previous ones
        analysis # This one
    )
}



analyses %>%
    glimpse()

# Before we write the raw "analyses" table, I want to add info from the database. This is going to make things easier in the report.
analyses %>%
    left_join(kegg_data %>% distinct(class, group, pathway), by = "pathway") %>%
    write_tsv(paste0(output_path, "/kegg_pathway_enrichment_analysis.tsv"))

analyses %>%
    select(sample, pathway, `p.adjust`) %>%
    left_join(kegg_data %>% distinct(class, group, pathway), by = "pathway") %>%
    select(sample, class, group, pathway, p_adj = `p.adjust`) %>%
    pivot_wider(names_from = sample, values_from = p_adj) %>%
    arrange(class, group, pathway) %>%
    write_tsv(paste0(output_path, "/summary.tsv"))
