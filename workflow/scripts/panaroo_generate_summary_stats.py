#!/usr/bin/env python

__author__ = "Carl M. Kobel"

import sys
import numpy as np 

input_file = sys.argv[1]

# This file reads the panaroo gene_presence_absence.Rtab file and outputs a roary/panaroo compatible summary statistics file. The reason why this is necessary is that panaroo does not output the summary_statistics.tsv file when the core genome is empty, but comparem2 needs that.
# This version only uses numpy in order to not impose any more dependencies into the docker image currently.

matrix = np.genfromtxt(fname=input_file, delimiter="\t", skip_header=1, filling_values=0)  # change filling_values as req'd to fill in missing values
#print(matrix)
#print("//")

n = len(matrix[0])-1 # The first row, minus the gene column.
m_pan = len(matrix) # number of genes (pan)

#print("n", n)
#print("m", m)
matrix_summed = matrix.sum(axis = 1)


m_core = sum(matrix_summed == n) # number of core genes

print(f"partition\tcount")
print(f"Core genes\t{m_core}")
print(f"Total genes\t{m_pan}")


