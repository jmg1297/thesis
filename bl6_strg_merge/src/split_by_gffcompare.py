'''
Script to create GTF files for each class created by gffcompare
'''

import sys
import os
import re
import argparse
import subprocess
from math import *

def get_gtf_line_dict(query_gtf):
    '''
    Create a dictionary linking transcript IDs to ordered lists of GTF lines
    '''
    print("Creating dictionary of query GTF lines")
    sys.stdout.flush()

    num_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'" % query_gtf, shell = True))
    percent = 0.0
    print('['+' '*50+']'),
    print('\r'),

    gtf_line_dict = {}
    with open(query_gtf, 'r') as gtf:
        for i,line in enumerate(gtf):

            percent = 100.0*(float(i)/num_lines)
            num_bars = int(floor(percent/2))
            print('\r'),
            print('['+'|'*num_bars+' '*(50 - num_bars) + ']'),
            sys.stdout.flush()

            if line[0] == "#":
                continue
            line_list = re.split("\s+", line.strip())
            transcript_id = line_list[11][1:-2]
            try:
                gtf_line_dict[transcript_id].append(line)
            except KeyError:
                gtf_line_dict[transcript_id] = [line]

    print('\r'),
    print('['+'|'*50 + ']'),

    print("\nDone")
    sys.stdout.flush()
    return gtf_line_dict

def split_gtf(gffcompare_tmap, query_gtf, out_dir):
    '''
    We will not create one file per class, but will group together classes of
    similar meaning/interest, as follows:
    '''
    file_to_class = {
        "annotated":["=", "c", "j", "e", "o"],
        "intronic":["i"],
        "novel_non-coding":["u", "x", "r"],
        "other_classes":["p", "s"]
    }

    class_to_file = {}
    for f in file_to_class.keys():
        for c in file_to_class[f]:
            class_to_file[c] = f

    # First we need to create a dictionary linking stringtie IDs to the relevant
    # lines from the query gtf, ready to be written to a file (and therefore
    # stored in order).
    gtf_line_dict = get_gtf_line_dict(query_gtf)

    outfile_dict = {f:open(os.path.join(out_dir, f+".gtf"), 'wa') for f in file_to_class.keys()}

    print("Reading through tmap file")
    sys.stdout.flush()

    num_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'" % gffcompare_tmap, shell = True)) - 1
    percent = 0.0

    with open(gffcompare_tmap, 'r') as tmap:
        # Skip the header line
        tmap.readline()
        for i,line in enumerate(tmap):

            percent = 100.0*(float(i)/num_lines)
            num_bars = int(floor(percent/2))
            print('\r'),
            print('['+'|'*num_bars+' '*(50 - num_bars) + ']'),
            sys.stdout.flush()

            line_list = re.split("\s+", line.strip())
            class_code = line_list[2]
            transcript_id = line_list[4]
            outfile_dict[class_to_file[class_code]].write("".join(gtf_line_dict[transcript_id]))

    print('\r'),
    print('['+'|'*50 + ']'),

    print("\nDone")
    sys.stdout.flush()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gffcompare_tmap")
    parser.add_argument("query_gtf")
    parser.add_argument("out_dir")
    args = parser.parse_args()

    split_gtf(args.gffcompare_tmap, args.query_gtf, args.out_dir)
