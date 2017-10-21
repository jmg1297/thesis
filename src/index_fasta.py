'''
Index a FASTA file and create an index file with following format:
first_line:last_line:name
'''

import sys


input_fa = sys.argv[1]
output_fname = input_fa + ".idx"

with open(input_fa, 'r') as f, open(output_fname, 'wa') as out:
    name = f.readline().strip().split()[0][1:]
    start = 1
    i = 1
    for line in f:
        i += 1
        if line[0] == ">":
            end = i - 1
            out.write(":".join(map(str, [start, end, name])) + "\n")
            name = line.strip().split()[0][1:]
            start = i

    end = i
    out.write(":".join(map(str, [start, end, name])) + "\n")
