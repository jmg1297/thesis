'''
Assign random transcription strands to a set of retrocopies, based on the
proportions of strand combinations found in expressed retrocopies, and write the
results to a BED file.
'''

import json
from numpy.random import choice as randchoice

def get_strand_props(strand_combo_counts_json):
    with open(strand_combo_counts_json, 'r') as f:
        strand_combo_counts = json.load(f)

    strands = ["+", "-"]
    strand_combo_props = {
        s:{
            t:{
                u:0.0 for u in strands
            } for t in strands
        } for s in strands
    }

    for ps, d in strand_combo_counts.iteritems():
        for rs, e in d.iteritems():
            total = sum(e.values())
            for ts, c in e.iteritems():
                strand_combo_props[ps][rs][ts] = float(c)/total
    print(strand_combo_props)
    return strand_combo_props

def assign_strands(mbed_fname, strand_combo_props, out_fname):
    strands = ["+", "-"]
    with open(mbed_fname, 'r') as f, open(out_fname, 'wa') as out:
        for line in f:
            line_list = line.strip().split()
            rc_strand = line_list[5]
            parent_strand = line_list[11]
            if parent_strand != ".":
                prob_dict = strand_combo_props[parent_strand][rc_strand]
                p = [prob_dict[s] for s in strands]
                transcript_strand = randchoice(strands, p=p)
                out_line = "\t".join(line_list[0:5] + [transcript_strand]) + "\n"
                out.write(out_line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mbed_fname", help="mBED file of retrocopies and their parents")
    parser.add_argument(
        "strand_combo_counts_json",
        help="JSON file of strand combinations and the count of each"
    )
    parser.add_argument("out_fname")
    args = parser.parse_args()

    strand_combo_props = get_strand_props(args.strand_combo_counts_json)

    assign_strands(args.mbed_fname, strand_combo_props, args.out_fname)
