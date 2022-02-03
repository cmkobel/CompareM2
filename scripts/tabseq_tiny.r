

read_fasta = function(file, basename_only = T, skip = 0) {

    # For spontaneous debugging:
    # file = "~/assemblycomparator2/tests/E._faecium_plasmids/VB3240.fna"

    if (basename_only) {
        file_presentable = basename(file)
    } else {
        file_presentable = file
    }

    # Open the file
    #cat(paste("reading", file, "as fasta", "\n"))
    message(paste("reading", file, "as fasta", "\n"))

    file_open = scan(file, what = "character", sep = "\n", quiet = T)


    rv = dplyr::tibble(raw = file_open[(skip+1):length(file_open)]) %>%

        # detect record headers
        dplyr::mutate(header = dplyr::if_else(stringr::str_sub(raw, 1, 1) == ">", T, F)) %>%

        # enumerate the record headers and fill the downwards through the sequence lines
        dplyr::mutate(header_num = ifelse(header, seq.int(nrow(.)), NA)) %>% # needs no sorting
        tidyr::fill(header_num, .direction = "down") %>%

        # Collect the lines for each record
        dplyr::group_by(header_num) %>%
        dplyr::summarize(part = stringr::str_sub(raw[1], 2),
                         sequence = paste(raw[-1], collapse = ""), .groups = "drop") %>%

        dplyr::mutate(comment = paste("record", dplyr::row_number(header_num))) %>%  # reset header nums

        dplyr::transmute(sample = file_presentable, part, comment, sequence)

    #cat(paste("parsed", (rv %>% dim)[1], "records", "\n"))
    message(paste("parsed", (rv %>% dim)[1], "records", "\n"))


    rv




}















write_fasta = function(x, file, format = "fasta", record_format = "%part", verbose = F) {

    if (format == "fasta") {

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

# TODO: The complement function should really be called iupac_complement









reverse_complement = function(string, reverse = TRUE, type = "nucleotide") {

    if (type == "nucleotide") {
        # Spontaneous debugging
        #string = "atgtcnNe"

        # Define a mapping 1:1
        mapping = c(
            # lower case
            "a" = "t", "t" = "a", "u" = "a", "c" = "g", "g" = "c", # lower case nucleotides
            "y" = "r", "r" = "y", "s" = "s", "w" = "w", "k" = "m", "m" = "k", # Two letter codes
            "b" = "v", "v" = "b", "d" = "h", "h" = "d", # Three letter codes
            "-" = "-", "n" = "n", # Miscellaneous

            # UPPER CASE
            "A" = "T", "T" = "A", "U" = "A", "C" = "G", "G" = "C", # LOWER CASE NUCLEOTIDES
            "Y" = "R", "R" = "Y", "S" = "S", "W" = "W", "K" = "M", "M" = "K", # TWO LETTER CODES
            "B" = "V", "V" = "B", "D" = "H", "H" = "D", # THREE LETTER CODES
            "-" = "-", "N" = "N") # MISCELLANEOUS


        #splitted = unlist(strsplit(string, ""))
        splitted = string %>% strsplit("") %>% unlist()
        rv = mapping[splitted] #%>% paste(collapse = "")

        NAs = sum(is.na(rv))
        if (NAs > 0) {
            warning(paste(NAs, "occurrences of unsupported charactors were complemented to a question mark:", paste(unique(splitted[is.na(rv)]), collapse = " ")))

            rv[is.na(rv)] = "?"
        }
    } else {
        stop("Only type \"nucleotide\" is supported for now.")
    }


    if (reverse) {
        rv = rev(rv)
    }

    paste(rv, collapse = "")
}

#"atugcyrswkmbdhvnatgatgatgatgATUGCYRSWKMBDHVNATGATGATGATG" %>% reverse_complement() %>% reverse_complement()
# attgcyrswkmbdhvnatgatgatgatgATTGCYRSWKMBDHVNATGATGATGATG










read_gff = function(file, parse_attributes = TRUE) {

    col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes") # Source: https://m.ensembl.org/info/website/upload/gff3.html

    # For debugging
    #file = "~/assemblycomparator2/tests/E._faecium_plasmids/output_asscom2/samples/VB3240/prokka/VB3240.gff"

    # Find the line number where the fasta file starts
    file_open = scan(file, what = "character", sep = "\n", quiet = T)
    annotation_fasta_split_line_number = str_detect(file_open, "^##FASTA") %>% which()
    # If there is a fasta part in the gff file, read it with tabseq.
    if(annotation_fasta_split_line_number %>% length() > 0) {
        write(paste("splitting the file between annotation and fasta at line number:", annotation_fasta_split_line_number), stderr())

        fasta = read_fasta(file, skip = annotation_fasta_split_line_number - 1) # Be informed that when you use the skip argument, it counts exclusive of comment lines. Fortunately, no comment argument has been implemented in tabseq::read_fasta().
    } else {
        warning("The gff file doesn't contain a fasta part. The fasta item in the return value will be `NA`")
        fasta = NA
    }

    annotation = read_tsv(file_open[1:(annotation_fasta_split_line_number - 1)], col_names = col_names, comment = "##") %>%
        # View()
        identity()

    list(annotation = annotation, fasta = fasta)
}



GC_content_nonvectorized = function(string, position = 0) {
    #string = "agcccatgtgaccagc"
    #string = "AbbCdd"

    # Uppercase as a means of homogenization
    string = string %>% toupper()

    # Each item in `splitted` is now a single character
    splitted = string %>%
        strsplit("") %>%
        unlist()


    if (position != 0 & position >= 1 & position <= 3) {
        #write(paste0("Info: Calculating GC", position, "."), stderr())
        splitted = splitted[(1:length(splitted)-1)%%3 == (position -1)]
    } else if (position == 0) {
        #write(paste0("Info: Calculating GC in all positions."), stderr())
    } else {
        stop("Invalid position. Please choose 0 (all positions) or 1 through 3 for GC1 through GC3, respectively.")
    }


    Gs = (splitted == "G") %>% sum()
    Cs = (splitted == "C") %>% sum()

    GCs = Gs + Cs

    GCs/length(splitted)

}













GC_content = function(string, position = 0) {
    message("calling vectorized edition")

    # Simply call the non vectorized function above, and pass the arguments.
    rv = Vectorize(GC_content_nonvectorized, vectorize.args = "string")(string, position)
    names(rv) = NULL
    rv
}



# TODO: Consider implementing a View function that strips or simplifies the sequence column.











ts_neutralize = function(input) {
    message("version x24")
    input %>%
        mutate(sequence = paste(str_sub(sequence, 1, 10), "(..)"))
}









tsView = function(input) {
    message("version x24")
    input %>%
        mutate(sequence = paste(str_sub(sequence, 1, 10), "(..)")) %>%
        View(title = paste("tabseq View"))
}




