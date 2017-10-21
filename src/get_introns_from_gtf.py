'''
Given a GTF file, get the introns of each transcript and write them to a BED
file.
'''

import re
import subprocess
import tempfile as tf

def get_attr_dict(attr_str):
    return dict(
        map(
            lambda y: re.split(' *"', y)[0:-1],
            map(
                lambda x: x.strip(),
                attr_str.split(";")[0:-1]
            )
        )
    )

def get_introns(gtf_fname, output_fname):
    output_lines = []
    exon_dict = {}

    with open(gtf_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split("\t")

            if line_list[2] not in ["transcript", "exon"]:
                continue

            attr_str = line_list[-1]
            attr_dict = get_attr_dict(attr_str)

            if line_list[2] == "transcript":
                exon_dict[attr_dict["transcript_id"]] = {
                    "chromosome":line_list[0],
                    "exon_bounds":[],
                    "strand":line_list[6]}

            elif line_list[2] == "exon":
                exon_dict[attr_dict["transcript_id"]]["exon_bounds"] \
                    += [int(line_list[3]), int(line_list[4])]

    for transcript, info in exon_dict.iteritems():
        exon_bounds = sorted(info["exon_bounds"])
        introns = [exon_bounds[1:-1][i:i+2] for i in xrange(0, len(exon_bounds[1:-1]), 2)]
        for i,intron in enumerate(introns):
            line = "\t".join(
                map(
                    str, ["chr"+info["chromosome"],
                            intron[0],
                            intron[1],
                            transcript+":intron{}".format(i+1),
                            0,
                            info["strand"]]))
            output_lines.append(line)

    with tf.NamedTemporaryFile() as tmp:
        tmp.write("\n".join(output_lines) + "\n")
        tmp.flush()
        subprocess.call(
            "sort -k1,1 -k2,2n {} > {}".format(tmp.name, output_fname),
            shell=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    get_introns(args.gtf_fname, args.output_fname)
