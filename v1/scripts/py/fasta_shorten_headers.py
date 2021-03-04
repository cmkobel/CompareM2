#!/usr/bin/env python3

# This program shortens the headers of a fasta file
# Pipe fasta file with STDIN, outputs with STDOUT
# It also serializes the lines

import sys



LENGTH = 30

i = 0
for line in sys.stdin:
    if line[0] == '>':
        i += 1
        new_header = str(line[1:(LENGTH-3)]).strip()
        print(f'>{i:03}{new_header}', end = '\n') 
    else:
        print(line, end = '') # Aleady contains a newline in the end.

