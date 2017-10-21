'''
Given a JSON file summarising alignments, create two CSV files containing
percentages and raw numbers for each sample.
'''

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("json_fname")
parser.add_argument("pc_out_fname")
parser.add_argument("num_out_fname")
args = parser.parse_args()

with open(args.json_fname, 'r') as f:
    data = json.load(f)

headers = ["sample", "uniq", "multimapped", "rRNA", "unmapped"]

with open(args.pc_out_fname, 'wa') as out:
    out.write(",".join(headers) + "\n")
    for sample in sorted(data["read_percents"].keys()):
        d = data["read_percents"][sample]
        l = [sample] + [str(d[h]) for h in headers[1:]]
        out_line = ",".join(l) + "\n"
        out.write(out_line)

with open(args.num_out_fname, 'wa') as out:
    out.write(",".join(headers) + "\n")
    for sample in sorted(data["read_nums"].keys()):
        d = data["read_nums"][sample]
        l = [sample] + [str(d[h]) for h in headers[1:]]
        out_line = ",".join(l) + "\n"
        out.write(out_line)
