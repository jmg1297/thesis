'''
Filter the results from bedtools intersect -wao between a set of exons and
the blocks of retrocopies based on reciprocal overlap as a proportion of total
exon/block length. Output a file in the bedtools intersect output format.
'''

import sys
import tempfile as tf
import subprocess

def get_transcript_dict(exon_fname):
    transcript_dict = {}
    with open(exon_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = ".".join(line_list[3].split(".")[0:3])
            start = int(line_list[1])
            end = int(line_list[2])
            length = end - start
            if transcript not in transcript_dict:
                transcript_dict[transcript] = {
                    "chromosome":line_list[0],
                    "start":start,
                    "end":end,
                    "strand":line_list[5],
                    "total_exon_length":length
                }
            else:
                transcript_dict[transcript]["start"] = min(start, transcript_dict[transcript]["start"])
                transcript_dict[transcript]["end"] = max(end, transcript_dict[transcript]["end"])
                transcript_dict[transcript]["total_exon_length"] += length
    return transcript_dict

def get_block_dict(block_fname):
    block_dict = {}
    with open(block_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            feature = line_list[3]
            start = int(line_list[1])
            end = int(line_list[2])
            if feature not in block_dict:
                block_dict[feature] = {
                    "chromosome":line_list[0],
                    "start":start,
                    "end":end,
                    "strand":line_list[5]
                }
            else:
                block_dict[feature]["start"] = min(start, block_dict[feature]["start"])
                block_dict[feature]["end"] = max(end, block_dict[feature]["end"])
    return block_dict

def get_intersect_dict(intersect_fname):
    column_names = [
        "exon_chromosome", "exon_start", "exon_end",
        "exon_name", "exon_score", "exon_strand",
        "block_chromosome", "block_start", "block_end",
        "block_name", "block_score", "block_strand",
        "overlap"
    ]
    intersect_dict = {}
    with open(intersect_fname, 'r') as f:
        for line in f:
            line_dict = {k:v for k,v in zip(column_names, line.strip().split())}

            if line_dict["block_start"] == "-1":
                continue

            transcript_name = ".".join(line_dict["exon_name"].split(".")[0:-1])
            element_name = line_dict["block_name"]


            if transcript_name in intersect_dict:
                try:
                    intersect_dict[transcript_name][element_name] \
                        += int(line_dict["overlap"])
                except KeyError:
                    intersect_dict[transcript_name][element_name] \
                        = int(line_dict["overlap"])
            else:
                intersect_dict[transcript_name] \
                    = {element_name:int(line_dict["overlap"])}

    return intersect_dict

def get_output_line(transcript, feature, transcript_dict, block_dict):
    transcript_list = [transcript_dict[transcript][k] for k in ["chromosome", "start", "end", "strand"]]
    transcript_list.insert(-1, transcript)
    transcript_list.insert(-1, 0)
    transcript_line = "\t".join(map(str, transcript_list))

    feature_list = [block_dict[feature][k] for k in ["chromosome", "start", "end", "strand"]]
    feature_list.insert(-1, feature)
    feature_list.insert(-1, 0)
    feature_line = "\t".join(map(str, feature_list))

    return (transcript_line, feature_line)

def filter_intersects(intersect_dict, transcript_dict, block_dict, overlap_cutoff):
    '''
    For each transcript, check whether (total overlap/total exon length) is
    greater than the overlap cutoff. If it is, check the same for the
    retrocopies it overlaps with. If this also exceeds the overlap_cutoff,
    put the pair in the output list.
    '''

    out_intersects = []

    for transcript, intersects in intersect_dict.iteritems():

        total_exon_length = transcript_dict[transcript]["total_exon_length"]

        for feature, overlap in intersects.iteritems():

            exon_prop = float(overlap)/total_exon_length

            if exon_prop >= overlap_cutoff:
                transcript_line, feature_line \
                    = get_output_line(transcript, feature, transcript_dict, block_dict)
                out_intersects.append([transcript_line, feature_line, str(overlap)])

    return out_intersects

def write_output_file(out_intersects, output_fname):
    with tf.NamedTemporaryFile() as tmp:
        for intersect in out_intersects:
            tmp.write("\t".join(intersect) + "\n")
        tmp.flush()

        subprocess.check_output(
            "sort -k1,1 -k2,2n {} > {}".format(tmp.name, output_fname),
            shell=True
        )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("intersect_fname")
    parser.add_argument("exon_fname")
    parser.add_argument("block_fname")
    parser.add_argument("overlap_cutoff", type=float)
    parser.add_argument("output_fname")
    args = parser.parse_args()

    transcript_dict = get_transcript_dict(args.exon_fname)
    block_dict = get_block_dict(args.block_fname)

    intersect_dict = get_intersect_dict(args.intersect_fname)

    out_intersects = filter_intersects(
        intersect_dict, transcript_dict, block_dict, args.overlap_cutoff)

    write_output_file(out_intersects, args.output_fname)
