# This file takes a fasta from stdin and prints out a tabseq to stdout with all of the lengths and GCs

message("loading tidyverse")
library(tidyverse)

development = F



args = commandArgs(trailingOnly=TRUE)
genome_fasta = args[1]

message("these are the args: ", args)
message("this is the genome: ", genome_fasta)


if (development) {

    genome_fasta = "../Downloads/AF31-27BH.fna" # This will be overwritten by the script.

    } else { # Do the following when NOT devoling

    message("Installing loading tabseq")
    # This is messy (mostly slow), but works.
    if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
    }

    devtools::install_github("cmkobel/tabseq")
}
# Then load the library in your script
library(tabseq)




genome_ro = tabseq::read_fasta(genome_fasta)


# Make a dummy table that contains a concatenation of everything
sum_of_parts = genome_ro |>
    group_by(sample) |>
    summarize(part = "ALL PARTS", comment = "concatenation of all contigs", sequence = paste(sequence, collapse = ""))



genome = genome_ro |>
    bind_rows(sum_of_parts) |>
    mutate(length = str_length(sequence),
           GC = tabseq::GC_content(sequence))



genome |>
    select(`#sample` = sample, part, length, GC) |>
    format_tsv() |>
    cat() # Write to stdout



