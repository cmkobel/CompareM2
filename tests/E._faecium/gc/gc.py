#!/usr/bin/env python

""" 
    This file takes all the genes from a fasta file and computes the 
    the mean GC3 content. The output is printed to STDOUT
"""

import sys, json

input_file = sys.argv[1]
try:
    cli_comment_1 = str(sys.argv[2])
    cli_comment_2 = str(sys.argv[3])
    
except Exception as e:
    print(e)
    cli_comment_1 = ""
    cli_comment2 = ""







def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)

dict = {}

#print('gene', 'posA', 'posB', 'comment', 'gc_content', 'cli_comment_1', 'cli_comment_2', sep = '\t')
print('header', 'length', 'GCs', 'length', 'GC3')

header = None

with open(input_file, 'r') as file:
    for line_raw in file:
        line = line_raw.strip()
        # if line[1] == "=":
        #     #print('skipping')
        #     del gene, right, pos, comment, posA, posB, dna, dna_gc, gc_content
        #     continue
        if line[0] == ">": # Header
            #eprint(line_s)

            # Dump last header
            if not header is None:
                print(header, length, GCs, length, GCs/float(length)/3, sep = '\t')

            
            header = f"\"{line[1:]}\""
            length = 0
            GCs = 0
            #eprint(header)

            
        else: # Subsequent DNA

            thirds = line[2::3].upper()
            
            GCs += sum([(i in ['G', 'C']) for i in thirds])
            length += len(thirds)


# Lastly, dump the last record
print(header, length, GCs, length, GCs/float(length)/3, sep = '\t')