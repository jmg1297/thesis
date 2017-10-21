'''
Script to get midpoints of mm10 retrocopies and their matches in CAST, and
create a circos connector data file
'''

import re

def get_ref_dict(retrocopy_bed):
    ref_dict = {}
    with open(retrocopy_bed, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            chrom = re.sub('chr', 'mm', line_list[0])
            midpoint = int(line_list[1]) + int(line_list[2])/2
            name = line_list[3]
            ref_dict[name] = {"chrom":chrom, "midpoint":midpoint, "matches":[]}

    return ref_dict

def add_matches(ref_dict, match_bed):
    with open(match_bed, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            names = line_list[3].split(",")
            midpoint = int(line_list[1]) + int(line_list[2])/2
            for name in names:
                ref_dict[name]["matches"].append(midpoint)

    return ref_dict

def write_output_lines(match_dict, output_fname):
    with open(output_fname, 'wa') as out:
        for name, info in match_dict.iteritems():
            for m in info["matches"]:
                line = " ".join(map(str, [info["chrom"], m, info["midpoint"]])) + "\n"
                out.write(line)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("retrocopy_bed")
    parser.add_argument("match_bed")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    ref_dict = get_ref_dict(args.retrocopy_bed)

    match_dict = add_matches(ref_dict, args.match_bed)

    write_output_lines(match_dict, args.output_fname)
