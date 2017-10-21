'''
Given a .annotated.gtf file from gffcompare, extract the lines matching a given
class code.
'''

import re

def get_attr_dict(attr_str):
    return dict(
        map(
            lambda y: re.split(' *"', y)[0:-1],
            map(
                lambda x: x.strip(),
                attr_str.split(";")[0:-1])))

def get_class_code(annotated_gtf, class_code, output_fname):
    keep_transcripts = []

    with open(annotated_gtf, 'r') as f, open(output_fname, 'wa') as out:
        for line in f:
            line_list = line.strip().split("\t")
            attr_dict = get_attr_dict(line_list[-1])

            if line[0:3] != "chr":
                line = "chr" + line

            if line_list[2] == "transcript":
                if attr_dict["class_code"] == class_code:
                    out.write(line)
                    keep_transcripts.append(attr_dict["transcript_id"])
            else:
                if attr_dict["transcript_id"] in keep_transcripts:
                    out.write(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("annotated_gtf")
    parser.add_argument("class_code")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    get_class_code(args.annotated_gtf, args.class_code, args.output_fname)
