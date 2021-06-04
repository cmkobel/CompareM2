library(tabseq)
library(tidyverse)



# The idea is to make a plot that shows the GC for each gene along the contigs.




# ~~~ Helper functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tview = function(table) {
    table |>
        rowwise() |>
        mutate(sequence = paste0(str_sub(sequence, 1, 7),
                                " ..(", max(0, str_length(sequence)-14), " bases).. ",
                                str_sub(sequence, -7))) |>
        View()
}

# A helper function that calculates the GC content in all positions
# Position 0 means that GC should be measured in all positions
# Positions 1 through 3 means that GC should be measured for GC1 through GC3, respectively.
GC_content = function(string, position = 0) {

    #string = "agcccatgtgaccagc"
    #string = "AbbCdd"

    string = string |> toupper()
    splitted = string |>
        strsplit("") |>
        unlist()


    if (position != 0 & position >= 1 & position <= 3) {
        #write(paste0("Info: Calculating GC", position, "."), stderr())
        splitted = splitted[(1:length(splitted)-1)%%3 == (position -1)]
        #write(paste0(splitted, collapse = ""), stderr())
    } else if (position == 0) {
        #write(paste0("Info: Calculating GC in all positions."), stderr())
    } else {
        stop("Invalid position. Please choose 0 (all positions) or 1 through 3 for GC1 through GC3, respectively.")
    }

    #GCs = (splitted == "G" | splitted == "C") |> sum()
    Gs = (splitted == "G") |> sum()
    Cs = (splitted == "C") |> sum()

    GCs = Gs + Cs

    GCs/length(splitted)

}

GC_content("AACTTG", position = 3)

assembly = tabseq::read_fasta("~/assemblycomparator2/tests/E._faecium_plasmids/VB3240.fna") |>
    #rowwise() %>%
    mutate(part = str_split(part, " ") |>
               map(1) |>
               unlist())

assembly |>
    tview()

assembly = assembly[1,]

annotation = read_tsv("~/assemblycomparator2/tests/E._faecium_plasmids/output_asscom2/samples/VB3240/prokka/VB3240.gff", comment = "##",
                      col_names = unlist(strsplit("seqid source type start end score strand phase attributes", split = " "))) |>
    rowwise() |>
    mutate(mid = mean(c(start, end), na.rm = T)) |>
    ungroup() |>

    mutate(biological_start = case_when(strand == "+" ~ start,
                                        strand == "-" ~ end)) |>
    drop_na(source)




        # Mark the fasta part
        mutate(fasta_record = case_when(str_detect(seqid, "^>") ~ TRUE)) |>
        fill(fasta_record)

    fasta_part = gff_df |>
        filter(is.na(fasta_record)) |>
        select(fasta_record) |>

        View()


        # if (parse_attributes) {
        mutate(attributes = str_split(attributes, ";")) |>
        unnest(attributes) |>View()
        separate(attributes, into = c("name", "value"), sep = "=") |>
        pivot_wider(seqid, names_from = "name", values_from = "value") |>

        View()

}


final = annotation |>
    left_join(assembly, by = c("seqid" = "part")) |>

    mutate(sequence = str_sub(sequence, start, end)) |>
    rowwise() |>
    mutate(sequence = case_when(strand == "+" ~ sequence,
                                strand == "-" ~ tabseq::reverse_complement(sequence))) |>

    mutate(GC3 = GC_content(sequence, position = 3)) |>
    #mutate(sequence = str_sub(sequence, start, end)) |>

    #mutate(modthree = (start - end) %%3) |>
    #tview()
    identity()

library(viridisLite)

mean_GC_value = final$GC3 |> mean(na.rm = T)
median_GC_value = final$GC3 |> median(na.rm = T)


final |> ggplot(aes(x = start, y = GC3, color = GC3)) +
    geom_hline(yintercept = mean_GC_value, linetype = "dashed", alpha = 0.5, size = 0.5) +


    #geom_smooth(span = .2, color = "red")+
    geom_smooth(span = .1, alpha = 0.2)+


    geom_segment(aes(x = start, xend = end, y = GC3, yend = GC3), size = 2) +
    geom_line(aes(x = biological_start), alpha = 0.4) +


    #scale_colour_viridis_b(direction = -1)+
    #scale_colour_gradient2(midpoint = mean_GC_value, low = rgb(.7, 0, 0), mid = rgb(0, 0.7, 0), high = rgb(0, 0, 0.7)) +
    scale_colour_gradient2(midpoint = mean_GC_value, low = rgb(.7, 0, 0), mid = rgb(0.7, 0.7, 0), high = rgb(0, 0.6, 0)) +

    scale_x_continuous(minor_breaks = seq(1, 10000000, 10000)) +

    #facet_grid(.~seqid, scales = "free_x") +
    theme_bw()

#ggsave("bigGC3.png", height = 10, width = 16)

ggsave("smallGC3.png", height = 2.5, width = 9)





GC_content("gggaaaaaa", position = 0)


