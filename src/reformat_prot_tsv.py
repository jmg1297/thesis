'''
Convert raw proteome TSVs to ones with the "Description" field split into its
constituent parts, and missing fields replaced with "NULL".
'''

import re
import sys

def quote(s):
    return "\"{}\"".format(s)

def reformat_tsv(input_tsv_fname, output_fname):
    out_lines = []
    with open(input_tsv_fname, 'r') as f:
        headers = [
            re.sub("\s+", "_", h) for h in f.readline().strip().split("\t")]

        raw_headers = headers[0:]

        metadata_keys = ["OS", "GN", "PE", "SV"]
        description_idx = headers.index("Description")
        for nh in reversed(metadata_keys):
            headers.insert(description_idx+1, nh)

        out_lines.append("\t".join(headers))

        for line in f:
            line_list = line.strip().split("\t")

            line_list = ["NULL" if x=="" else x for x in line_list]

            line_dict = {h:l for h,l in zip(raw_headers, line_list)}

            full_description = line_dict["Description"]
            description = quote(" ".join(full_description.split("=")[0].split()[0:-1]))
            metadata_pattern = re.compile(r'([A-Z]{2}=.*)+')
            metadata = re.search(metadata_pattern, full_description).group()
            rev_items_pattern = re.compile(r'.+?=[A-Z]{2} ?')

            items = dict([
                tuple(
                    map(
                        lambda x:re.sub("\s+", "_", x.strip()),
                        i.split("=")))
                for i in reversed(
                    [s[::-1] for s in re.findall(
                        rev_items_pattern, metadata[::-1])])])

            out_line = "\t".join(
                line_list[0:description_idx] \
                + [description] \
                + [items.get(k, "NULL") for k in metadata_keys] \
                + line_list[description_idx+1:])

            out_lines.append(out_line)

    with open(output_fname, 'wa') as out:
        out.write("\n".join(out_lines))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    reformat_tsv(args.input_tsv_fname, args.output_fname)
