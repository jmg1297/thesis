'''
Filter tabular output of blastn to leave only those where the query
chromosome matches the subject chromosome, and the match length as a proportion
of the query length is above a certain threshold.
'''

import logging
from collections import OrderedDict

logging.basicConfig(level=logging.INFO)

FIELDS=OrderedDict([
    ("query_name",str),
    ("subject_name",str),
    ("percent_identity",float),
    ("aln_length",int),
    ("mismatches",int),
    ("gap_opens",int),
    ("query_start",int),
    ("query_end",int),
    ("subject_start",int),
    ("subject_end",int),
    ("evalue",float),
    ("bit_score",float),
    ("line",str)
])

def line_to_dict(line):
    line_list = line.strip().split()
    return {f:t(l) for (f,t),l in zip(FIELDS.iteritems(), line_list + [line])}

def get_blastn_out_dict(blast_out_fname):
    logging.info("Loading blastn output")
    blastn_out_dict = {}
    with open(blast_out_fname, 'r') as f:
        for i,line in enumerate(f):
            if (i+1)%1000 == 0:
                logging.info("{} lines parsed".format(i+1))

            line_list = line.strip().split()
            name = line_list[0]
            try:
                blastn_out_dict[name].append(line_to_dict(line))
            except KeyError:
                blastn_out_dict[name] = [line_to_dict(line)]
    logging.info("All lines parsed, returning dict")
    return blastn_out_dict

def filter_blastn(blastn_out_dict, query_bed_fname, min_length_prop):
    logging.info("Filtering blastn output")

    filtered = []

    with open(query_bed_fname, 'r') as f:
        for i,line in enumerate(f):
            if (i+1)%1000 == 0:
                logging.info("{} query lines parsed".format(i+1))

            line_list = line.strip().split()
            chrom = line_list[0]
            length = int(line_list[2]) - int(line_list[1])
            name = line_list[3]

            try:
                matches = blastn_out_dict[name]
            except KeyError:
                continue

            for match in matches:
                if chrom != match["subject_name"]:
                    continue
                else:
                    length_prop = float(match["aln_length"])/length
                    if length_prop < min_length_prop:
                        continue
                    else:
                        filtered.append(match["line"])

    logging.info("Finished filtering")
    return filtered

def write_output(filtered, output_fname):
    logging.info("Writing output")
    with open(output_fname, 'wa') as out:
        out.write("".join(filtered))
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_out_fname")
    parser.add_argument("query_bed_fname")
    parser.add_argument("min_length_prop", type=float)
    parser.add_argument("output_fname")
    args = parser.parse_args()

    blastn_out_dict = get_blastn_out_dict(args.blast_out_fname)

    filtered = filter_blastn(
        blastn_out_dict, args.query_bed_fname, args.min_length_prop)

    write_output(filtered, args.output_fname)
