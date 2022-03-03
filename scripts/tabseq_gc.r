# This file takes a fasta from stdin and prints out a tabseq to stdout with all of the lengths and GCs



message("loading tidyverse")
library(tidyverse)

development = F



args = commandArgs(trailingOnly=TRUE)
tiny_tabseq = args[1]
genome_fasta = args[2]

message("these are the args: ", args)
message("this is tiny tabseq: ", tiny_tabseq)
message("this is the genome: ", genome_fasta)


message("sourcing tiny tabseq")
source(tiny_tabseq)


if (development) {

    genome_fasta = "../Downloads/AF31-27BH.fna" # This will be overwritten by the script.
    source("../assemblycomparator2/scripts/tabseq_tiny.r")

} else { # Do the following when NOT devoling
    a = 3 # equivalent to doing nothing
    #message("Installing loading tabseq")
    # This is messy (mostly slow), but works.
    # if (!requireNamespace("devtools", quietly = TRUE)) {
    #     install.packages("devtools", repos = "http://cran.us.r-project.org")
    # }

    #devtools::install_github("cmkobel/tabseq@main")
}
# Then load the library in your script
#library(tabseq)




genome_ro = read_fasta(genome_fasta)


# Make a dummy table that contains a concatenation of everything
sum_of_parts = genome_ro %>%
    group_by(sample) %>%
    summarize(part = "ALL PARTS", comment = "concatenation of all contigs", sequence = paste(sequence, collapse = ""))



genome = genome_ro %>%
    bind_rows(sum_of_parts) %>%
    mutate(length = str_length(sequence),
           GC = GC_content(sequence))




genome %>%
    select(`#sample` = sample, part, length, GC) %>%
    format_tsv() %>%
    cat()


