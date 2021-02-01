#!/usr/bin/env python3


import sys
import argparse

__author__ = "Carl Mathias Kobel"
__version__ = "0.1.1" # Not tested enough


def check_non_negative_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue



# argparser stuff
parser = argparse.ArgumentParser(prog = "tabseq2fasta.py", description='Convert .tabseq to .fasta from stdin to stdout')
parser.add_argument('--version', action='version', version=f"tabseq2fasta.py v{__version__}")

# Options to modify data
# TODO: In order to implement this, there should be a mutually exclusive group for each variable, and I feel like it would complicate things unnecessarily.
#parser.add_argument('--fill_sample', nargs='?', help='overwrite the sample data with this value', default = False)
#parser.add_argument('--fill_part', nargs='?', help='overwrite the part data with this value', default = False)
#parser.add_argument('--fill_comment', nargs='?', help='overwrite the comment data with this value', default = False)

# Options to remove data
group = parser.add_argument_group('optionally, omit data from being written to the fasta records')
group.add_argument("--clear_sample", help = "clear sample data", action = "store_true")
group.add_argument("--clear_part", help = "clear part data", action = "store_true")
group.add_argument("--clear_comment", help = "clear comment data", action = "store_true")
group.add_argument("--clear_all", help = "shortcut for the above clear-arguments", action = "store_true")

# Option to enumerate the records, and even specify the number of padded zeros
parser.add_argument("--enumerate", help = "prefix the sequences with numbers, optionally specify the number of padded zeros, defaults to 3 padded zeros", nargs = "?", type = check_non_negative_int, metavar = "number of zeros", default = False, const = 3)

# Option to wrap the fasta file lines
parser.add_argument("--wrap", help = "wrap the fasta output to a limited line length, defaults to 60 characters", nargs = "?", type = check_non_negative_int, metavar = "line length", default = False, const = 60)

# 
parser.add_argument("--separator", help = "choose a symbol to separate the fasta records, defaults to \"|\"", default = "|")



args = parser.parse_args()



def eprint(string, *args, **kwargs):
    print(string, *args, **kwargs, file = sys.stderr)


f"""
    This is the tabseq2fasta_cf.py script
    It ouputs a wrapped fasta file to output.
    The script might need testing?
"""





def wrap(sequence, width):
    # Wrap
    rv = ''
    for i in range(0, len(sequence), width):
        rv += sequence[i:(i+width)] + '\n'
    return rv.strip()

num = 0


for _i, line in enumerate(sys.stdin):
    if line[0] == '#': # Skip headers 
        continue

    sample, part, comment, sequence = line.split('\t')
    sequence = sequence.strip("\n")


    if args.clear_all:
        sample, part, comment = "", "", ""
    else:
        if args.clear_sample: sample = ""
        
        if args.clear_part:
            part = ""
        else: part = args.separator + part
        
        if args.clear_comment:
            comment = ""
        else: comment = args.separator + comment

    if not args.enumerate:
        #enum = f"{_i:01}{args.separator}"
        enum = ""
    else:
        enum = f"{str(_i).rjust(args.enumerate, '0')}{args.separator}"
    


    fasta_header = f">{enum}{sample}{part}{comment}" # "look at the tree" # A_NA_3328-130A
    #eprint(fasta_header[1:])
    print(fasta_header)


    try:
        if not args.wrap:
            print(sequence)
        else:
            print(wrap(sequence, args.wrap))
    except BrokenPipeError as e: # Gracefully handles shell handling, i.e the head command.
        pass

    num += 1





eprint('\nWrote', num, 'sequences to stdout.\n')