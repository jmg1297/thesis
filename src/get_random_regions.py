'''
Script to generate a BED file of randomly chosen regions, according to a length
distribution obtained from a reference BED file, as well as the ref dist across
chromosomes
'''

import sys
import random
from numpy.random import choice as npchoice

def jitter(length, jitter_size):
    return random.randint(length-jitter_size, length+jitter_size)

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
            chrom_sizes[chrom] = int(size)

    print("Done")
    sys.stdout.flush()

    return chrom_sizes

def get_random_regions(chrom_sizes, ref_fname, jitter_size, sample_size, out_fname):

    print("Getting reference length and chromosome distributions ... "),
    sys.stdout.flush()

    chrom_counts = {}
    ref_lengths = []
    with open(ref_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            try:
                chrom_counts[line_list[0]] += 1
            except KeyError:
                chrom_counts[line_list[0]] = 1
            ref_lengths.append(int(line_list[2]) - int(line_list[1]))

    ref_total = sum(chrom_counts.values())
    chrom_props = {
        chrom:float(count)/ref_total for chrom,count in chrom_counts.iteritems()
    }

    chroms = sorted(chrom_props.keys())
    probs = [chrom_props[c] for c in chroms]

    print("Done")

    print("Generating random chromosomes and lengths ... "),
    sys.stdout.flush()

    sample_chroms = npchoice(chroms, size=sample_size, p=probs)
    sample_lengths = map(
        lambda x: jitter(x, jitter_size),
        npchoice(ref_lengths, size=sample_size).tolist()
    )

    print("Done")
    print("Writing random regions to file ... "),
    sys.stdout.flush()

    with open(out_fname, 'wa') as out:
        for i, (chrom, length) in enumerate(zip(sample_chroms, sample_lengths)):
            chrom_size = chrom_sizes[chrom]
            start = random.randint(1, chrom_size-length)
            end = start + length
            name = "RAND"+str(i)
            score="1"
            strand="+"
            out.write(
                "\t".join([chrom, str(start), str(end), name, score, strand]) + "\n")

    print("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("chrom_sizes_fname")
    parser.add_argument("ref_fname")
    parser.add_argument("jitter_size", type=int)
    parser.add_argument("sample_size", type=int)
    parser.add_argument("out_fname")
    args = parser.parse_args()

    chrom_sizes = get_chrom_size_dict(args.chrom_sizes_fname)

    get_random_regions(
        chrom_sizes,
        args.ref_fname,
        args.jitter_size,
        args.sample_size,
        args.out_fname
    )
