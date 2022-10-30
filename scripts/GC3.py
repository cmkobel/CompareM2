#!/usr/bin/env python3

""" 
    This file takes all the genes from a fasta file and computes the 
    the mean GC3 content. The output is printed to STDOUT
"""

import sys, json


try:
    input_file = str(sys.argv[1])
except IndexError as e:
    print(f"IndexError: {e}.\nPlease provide a fasta file as argument to the {sys.argv[0]} script.")
    print(f"\
Usage: \n\
    {sys.argv[0]} path/to/codon_aligned_genes.fasta [sample_name]\n\
\n\
")
    exit()






try:
    sample_name = str(sys.argv[2])    
except IndexError as e:
    sample_name = "NA"







def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)




#print('#sample', 'header', 'length', 'n_GC3', 'GC3', sep = '\t')



def calculate_and_print(header, DNA):
    length = len(DNA)
    thirds = DNA[2::3].upper()
    GCs = sum([(i in ['G', 'C']) for i in thirds])

    print(sample_name, header, length, GCs, (GCs/float(length))*3, sep = '\t')

header = None

with open(input_file, 'r') as file:
    for line_raw in file:
        line = line_raw.strip()
        # if line[1] == "=":
        #     #print('skipping')
        #     del gene, right, pos, comment, posA, posB, dna, dna_gc, gc_content
        #     continue
        if line[:1] == ">": # Header

            # Dump last header
            if not header is None:
                calculate_and_print(header, DNA)

            # Save new header and reset DNA string
            header = line[1:].strip().replace("\t", "\\t") # Make sure that tabs in the header string doesn't break the output tsv file.
            DNA = ""

            
        else: # Subsequent DNA (possibly multiline)

            DNA += line

            


# Lastly, dump the last record
calculate_and_print(header, DNA)