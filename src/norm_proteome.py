'''
Given a TSV file of proteome results, apply a log2 transformation and median
normalisation. Procedure:
- for each abundance column, get the median
- get the mean of the medians (M)
- get the ratio of each median to M
- divide each abundance value by the appropriate ratio
- log2 each abundance
'''

import sys
import numpy as np

def get_median_ratios(proteome_fname, abundance_fields):

    values = {f:[] for f in abundance_fields}

    with open(proteome_fname, 'r') as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            line_dict = {h:v for h,v in zip(headers, line.strip().split("\t"))}
            for field in abundance_fields:
                try:
                    values[field].append(float(line_dict[field]))
                except ValueError:
                    continue

    medians = {f:np.median(values[f]) for f in abundance_fields}

    median_mean = np.mean(medians.values())

    median_ratios = {f:medians[f]/median_mean for f in abundance_fields}

    return median_ratios

def transform(value, median_ratio):
    if value == "NULL":
        return "NULL"
    else:
        value = float(value)

    if value == 0:
        return 0
    else:
        return np.log2(100.0*(value/median_ratio))

def normalise_data(proteome_fname, median_ratios, output_fname):
    output_lines = []
    with open(proteome_fname, 'r') as f:
        header_line = f.readline()
        output_lines.append(header_line.strip())
        headers = header_line.strip().split("\t")
        for line in f:
            line_dict = {h:v for h,v in zip(headers, line.strip().split("\t"))}
            for f in median_ratios.keys():
                line_dict[f] = transform(line_dict[f], median_ratios[f])
            out_line = "\t".join(map(str, [line_dict[h] for h in headers]))
            output_lines.append(out_line)

    with open(output_fname, 'wa') as out:
        out.write("\n".join(output_lines) + "\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("proteome_fname")
    parser.add_argument("-f", "--abundance_field", action="append")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    median_ratios = get_median_ratios(args.proteome_fname, args.abundance_field)

    normalise_data(args.proteome_fname, median_ratios, args.output_fname)
