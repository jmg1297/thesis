'''
Get ranges of possible TSS from output of bedtools cluster for clusters
with more than sample contributing to them.
'''

import re
import subprocess
import tempfile as tf
import numpy as np

def get_cluster_dict(cluster_fname):
    cluster_dict = {}
    with open(cluster_fname, 'r') as f:
        for line in f:
            cluster = re.split("\s+", line.strip())[-1]
            try:
                cluster_dict[cluster].append(line)
            except KeyError:
                cluster_dict[cluster] = [line]
    return cluster_dict

def filter_clusters(cluster_dict):
    filtered_clusters = {}
    for cluster,lines in cluster_dict.iteritems():
        samples = set(
            [re.split("\s+", l.strip())[3].split(":")[0] for l in lines]
        )
        if len(samples) > 1:
            filtered_clusters[cluster] = lines
    return filtered_clusters

def get_consensus_boundaries(starts, ends):
    labelled_starts = [(s, "s") for s in starts]
    labelled_ends = [(e, "e") for e in ends]
    labelled_boundaries = sorted(
        labelled_starts + labelled_ends, key=lambda x: x[0]
    )
    boundaries = [x[0] for x in labelled_boundaries]
    labels = [x[1] for x in labelled_boundaries]
    consensus_starts = [boundaries[0], boundaries[labels.index('e')-1]]
    consensus_ends = [
        boundaries[len(labels) - labels[::-1].index('s')],
        boundaries[-1]
    ]

    return consensus_starts, consensus_ends

def get_tss_line(cluster, lines):
    chrom = re.split("\s+", lines[0].strip())[0]
    starts = []
    ends = []
    tpms = []
    strands = []
    for line in lines:
        line_list = re.split("\s+", line.strip())
        starts.append(int(line_list[1]))
        ends.append(int(line_list[2]))
        tpms.append(float(line_list[4]))
        strands.append(line_list[5])

    if len(set(strands)) != 1:
        consensus_strand = "."
    else:
        consensus_strand = list(set(strands))[0]

    consensus_starts, consensus_ends = get_consensus_boundaries(starts, ends)
    if consensus_strand == "-":
        tss_range = consensus_ends
    else:
        tss_range = consensus_starts

    tss_line = "\t".join(
        [chrom] \
        + map(str, tss_range)
        + ["CLUS"+str(cluster), str(np.mean(tpms)), consensus_strand]
    ) + "\n"
    print(tss_line),
    return tss_line

def write_tss_range_bed(filtered_clusters, out_fname):
    with tf.NamedTemporaryFile() as tmp:
        for cluster,lines in filtered_clusters.iteritems():
            tss_line = get_tss_line(cluster, lines)
            tmp.write(tss_line)
            tmp.flush()

        subprocess.check_output(
            "sort -k1,1 -k2,2n %s > %s" % (tmp.name, out_fname),
            shell=True
        )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("cluster_fname")
    parser.add_argument("out_fname")
    args = parser.parse_args()

    cluster_dict = get_cluster_dict(args.cluster_fname)
    filtered_clusters = filter_clusters(cluster_dict)
    write_tss_range_bed(filtered_clusters, args.out_fname)
