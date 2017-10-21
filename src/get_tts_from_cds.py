'''
Use a BED file of CDS regions to get transcript TSSs
'''

import sys
import tempfile as tf
import subprocess

input_bed = sys.argv[1]
out_bed = sys.argv[2]

with open(input_bed, 'r') as f:
    tts_dict = {}
    # Assume all CDS regions for a transcript are together
    line_list = f.readline().strip().split()
    current_id = line_list[3]
    current_chrom = line_list[0]
    current_strand = line_list[5]
    bounds = [int(line_list[1]), int(line_list[2])]

    for line in f:
        line_list = line.strip().split()
        if line_list[3] == current_id:
            bounds += [int(line_list[1]), int(line_list[2])]
        else:
            if current_strand == "+":
                tts = max(bounds)
            elif current_strand == "-":
                tts = min(bounds)
                
            tts_dict[current_id] = {
                "tts":tts, "chrom":current_chrom, "strand":current_strand}

            current_id = line_list[3]
            current_chrom = line_list[0]
            current_strand = line_list[5]
            bounds = [int(line_list[1]), int(line_list[2])]


with tf.NamedTemporaryFile() as tmp:
    for transcript_id, info in tts_dict.iteritems():
        line = "\t".join(
            [
                info["chrom"],
                str(info["tts"]),
                str(info["tts"]),
                transcript_id,
                "0",
                info["strand"]
            ]
        ) + "\n"
        tmp.write(line)
    tmp.flush()
    subprocess.check_output(
        "sort -k1,1 -k2,2n {} > {}".format(tmp.name, out_bed),
        shell=True
    )
