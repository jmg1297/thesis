'''
Given an mBED file of transcripts with intersecting retrotransposon elements,
write out a file with one line per transcript with a summary of overlaps for
each.
'''

import tempfile as tf
import subprocess

def get_output_lines(intersect_dict):
    output_lines = []
    for transcript, info in intersect_dict.iteritems():
        overlaps = []
        for element_type, overlap in info["overlaps"].iteritems():
            overlaps.append("{},{}".format(element_type, overlap))
        line = "\t".join([info["line"], ";".join(overlaps)])
        output_lines.append(line)
    return output_lines

def summarise_intersects(mbed_fname, output_fname):
    intersect_dict = {}

    with open(mbed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = line_list[3]
            element_type = ":".join(line_list[9].split(":")[0:2])
            overlap = int(line_list[-1])
            if transcript not in intersect_dict:
                transcript_line = "\t".join(line_list[0:6])
                intersect_dict[transcript] = {"line":transcript_line, "overlaps":{element_type:overlap}}
            else:
                try:
                    intersect_dict[transcript]["overlaps"][element_type] += overlap
                except KeyError:
                    intersect_dict[transcript]["overlaps"][element_type] = overlap

    output_lines = get_output_lines(intersect_dict)

    with tf.NamedTemporaryFile() as tmp:
        for line in output_lines:
            tmp.write(line + "\n")
        tmp.flush()

        subprocess.check_output(
            "sort -k1,1 -k2,2n {} > {}".format(tmp.name, output_fname),
            shell=True
        )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mbed_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    summarise_intersects(args.mbed_fname, args.output_fname)
