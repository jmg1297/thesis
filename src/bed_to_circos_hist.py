'''
Use a BED file to create a data file that can be used by Circos to create a
histogram.
'''

import sys
import os
import re
import argparse
from math import *
import subprocess

def get_chrom_size_dict(chrom_sizes_fname):
    '''
    Process a file of chromosome sizes to produce a dictionary of chromosome
    sizes
    '''
    print("Obtaining chromosome sizes ... "),
    sys.stdout.flush()

    chrom_sizes = {}
    with open(chrom_sizes_fname, 'r') as f:
        for line in f:
            chrom, size = line.strip().split()
            if "chr" in chrom:
                chrom = chrom[3:]
            chrom_sizes[chrom] = int(size)

    print("Done")
    sys.stdout.flush()

    return chrom_sizes

def init_count_dict(chrom_sizes, bin_size):
    print("Setting up count dictionary ... "),
    sys.stdout.flush()

    count_dict = {chrom:{} for chrom in chrom_sizes.keys()}

    for chrom in count_dict.keys():
        size = chrom_sizes[chrom]
        num_bins = float(size)/bin_size
        if num_bins != int(num_bins):
            num_bins = int(ceil(num_bins))
        else:
            num_bins = int(num_bins)

        count_dict[chrom] = {i:{"start":0, "end":0, "value":0} for i in range(num_bins)}

        for i in range(num_bins):
            count_dict[chrom][i]["start"] = i*bin_size
            if (i+1)*bin_size < chrom_sizes[chrom]:
                count_dict[chrom][i]["end"] = (i+1)*bin_size
            else:
                count_dict[chrom][i]["end"] = chrom_sizes[chrom]

    print("Done")
    sys.stdout.flush()

    return count_dict

def write_circos(data_dict, chrom_prefix, out_fname):
    '''
    Use a dictionary containing chromosome, start, end, and histogram height
    data to write a data file that can be used by circos
    '''

    print("Writing circos file to {}".format(out_fname))
    sys.stdout.flush()

    with open(out_fname, 'wa') as f:
        for chrom in data_dict.keys():
            circos_chrom = chrom_prefix + chrom
            for i in sorted(data_dict[chrom].keys()):
                start = str(data_dict[chrom][i]["start"])
                end = str(data_dict[chrom][i]["end"])
                value = data_dict[chrom][i]["value"]
                f.write(
                    " ".join([circos_chrom, start, end, "%.2f" % value]) + "\n"
                )

def bed_to_circos_hist(chrom_sizes, bed_fname, chrom_prefix, bin_size, out_fname):
    '''
    Parse a BED file of regions and create a circos data file that can be used
    to create a circos histogram showing number of regions in each bin
    '''

    count_dict = init_count_dict(chrom_sizes, bin_size)

    num_lines = int(
        subprocess.check_output(
            "wc -l %s | awk '{print $1}'" % bed_fname,
            shell = True
        )
    )
    percent = 0
    print('['+' '*50+']'),
    print('\r'),

    with open(bed_fname, 'r') as f:
        for i,line in enumerate(f):
            line_list = re.split("\s+", line.strip())
            chrom = line_list[0]
            if "chr" in chrom:
                chrom = chrom[3:]
            start = int(line_list[1])
            end = int(line_list[2])

            if chrom not in chrom_sizes.keys():
                continue

            while True:
                start_bin = int(start/bin_size)
                # If the start and end of the region lie in the same bin, the
                # whole thing contributes to this value
                if end <= count_dict[chrom][start_bin]["end"]:
                    count_dict[chrom][start_bin]["value"] += 1
                    break
                # Otherwise, add the proportion of base pairs included in the start bin, then
                # update the start value to the start of the next bin and carry on
                else:
                    count_dict[chrom][start_bin]["value"] \
                        += float(count_dict[chrom][start_bin]["end"] - start)/(end-start)
                    start = count_dict[chrom][start_bin+1]["start"]

        percent = 100.0*(float(i)/num_lines)
        num_bars = int(floor(percent/2))
        print('\r'),
        print('['+'|'*num_bars+' '*(50-num_bars)+']'),
        sys.stdout.flush()

    print('\r'),
    print('['+'|'*50+']'),
    sys.stdout.flush()

    print("\nFinished going through file {}".format(bed_fname))
    sys.stdout.flush()

    print("Calculating relative values for circos plot ... "),
    sys.stdout.flush()

    '''
    circos_data = count_dict.copy()

    for chrom in circos_data.keys():
        for i in circos_data[chrom].keys():
            circos_data[chrom][i]["value"] \
                = 1000000*float(circos_data[chrom][i]["value"])/chrom_sizes[chrom]
    '''

    print("Done")
    print("Finished, writing circos data")
    sys.stdout.flush()

    write_circos(count_dict, chrom_prefix, out_fname)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("chrom_sizes_fname")
    parser.add_argument("bed_fname")
    parser.add_argument("chrom_prefix")
    parser.add_argument("bin_size", type=int)
    parser.add_argument("out_fname")
    args = parser.parse_args()

    chrom_sizes = get_chrom_size_dict(args.chrom_sizes_fname)
    bed_to_circos_hist(chrom_sizes, args.bed_fname, args.chrom_prefix, args.bin_size, args.out_fname)
