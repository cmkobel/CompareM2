#' Read a fasta file and convert it to tabseq
#' @export
#' @import dplyr stringr tidyr readr
#' @description Takes the path to a fasta file. Loads the fasta file, and converts it to a tabseq table.
#' @param file The path to a fasta-file
#' @param from_fasta A boolean designating whether the file path points to a fasta-file
#' @return A tibble containing the respective records and sequences from the fasta-file. Represented in the tabseq format.
#' @examples
#'
#' #my_sequences = read_tabseq("path/to/sequences.fa")
#'
read_tabseq = function(file, from_fasta = T) {


    if (from_fasta) {
        # Open the file
        cat(paste("reading", file, "as fasta", "\n"))
        file_open = scan(file, what = "character", sep = NULL, quiet = T)


        rv = dplyr::tibble(raw = file_open) %>%

            # detect record headers
            dplyr::mutate(header = if_else(stringr::str_sub(raw, 1, 1) == ">", T, F)) %>%

            # enumerate the record headers and fill the downwards through the sequence lines
            dplyr::mutate(header_num = ifelse(header, seq.int(nrow(.)), NA)) %>% # needs no sorting
            tidyr::fill(header_num, .direction = "down") %>%

            # Collect the lines for each record
            dplyr::group_by(header_num) %>%
            dplyr::summarize(part = stringr::str_sub(raw[1], 2),
                             sequence = paste(raw[-1], collapse = ""), .groups = "drop") %>%

            dplyr::mutate(comment = paste("record", dplyr::row_number(header_num))) %>%  # reset header nums

            dplyr::transmute(sample = file, part, comment, sequence)

        cat(paste("parsed", (rv %>% dim)[1], "records", "\n"))

        rv


    }

}

#' Write a tabseq table to disk in fasta-format.
#' @export
#' @import dplyr tidyr stringr
#' @description Takes a tabseq table, and saves it to disk in the fasta format.
#' @param x The tabseq table to save
#' @param file The path to a wanted fasta-file
#' @param to_fasta A boolean designating whether the file should be written in the fasta format.
#' @param record_format A string designating how to encode the records in the fasta file written. Defaults to "\%part", but may be a combination of all metadata columns including custom separators i.e. "\%sample|\%part|\%comment.
#' @param verbose For debugging
#' @return Saves the content to disk, and returns a tibble which represents what is being written to disk.
#' @examples
#'
#' #my_sequences = read_tabseq("path/to/sequences.fa")
#' #write_tabseq(my_sequences, "write/sequences/here/sequences.fa")
write_tabseq = function(x, file, record_format = "%part", to_fasta = T, verbose = F) {

    if (to_fasta) {


        cat(paste("writing to file ", file, "as fasta", "\n"))
        cat(paste("using the following record_format:", record_format, "\n"))


        # The only reason I'm placing the record definition before making the 'formatted' variable, is because I want to be able to verbosely read the contents of it.
        record = stringr::str_replace_all(record_format, "%sample", as.character(x$sample)) %>%
            stringr::str_replace_all("%part", as.character(x$part)) %>%
            stringr::str_replace_all("%comment", as.character(x$comment))

        if (verbose) {
            print(record)
        }


        # format as "fasta"
        formatted = x %>%
            dplyr::mutate(record_format = record_format) %>%
            dplyr::transmute(record = record, sequence = sequence)





        # Interleave
        interleave = as.vector(rbind(paste0(">", formatted$record), formatted$sequence))


        fileConn<-file(file)
        writeLines(interleave, fileConn)
        close(fileConn)

        #cat("\n")
        #cat(paste0(">", formatted$record, "\t", str_sub(formatted$sequence, 1, 10), "...", collapse = "\n"))
        formatted
    }

}




