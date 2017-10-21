'''
Given a list of transcript IDs, produces a methylation plot for each of the
corresponding regions. Each plot will show data from each of a given set of
methylation bedgraphs. The region can be extended in both directions by a given
amount. There will also be an extra plot showing the structure and direction of
the transcript.
'''

import sys
import os
import traceback as tb
import re
import linecache
import logging
import multiprocessing as mp
from Queue import Empty
from math import floor

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

logging.basicConfig(level=logging.INFO)

class GetMethProcess(mp.Process):
    def __init__(self, job_queue, results_queue):
        mp.Process.__init__(self)
        self.job_queue = job_queue
        self.results_queue = results_queue
        logging.info("Initialised {}".format(self.name))

    def run(self):
        logging.info("{} starting to retrieve jobs from queue".format(self.name))

        while True:
            job = self.job_queue.get()
            if job == "STOP":
                logging.info("Received STOP signal, stopping {}".format(self.name))
                return 0
            else:
                t, l, trans_info, add_left, add_right, bg, idx = job
                logging.info("Retrieved {}/{} pair, processing".format(t, l))
                cpgs = self.get_meth_data(trans_info, add_left, add_right, bg, idx)
                logging.info("Finished {}/{} pair, putting results on queue".format(t, l))
                self.results_queue.put([t, l, cpgs])

    def get_meth_data(self, trans_info, add_left, add_right, bedgraph_fname, idx_dict):

        chrom = trans_info["chromosome"]
        region_start = trans_info["start"] - add_left
        region_end = trans_info["end"] + add_right

        first_line, last_line = self.bisection_search(
            bedgraph_fname,
            idx_dict[chrom]["first"],
            idx_dict[chrom]["last"],
            region_start,
            region_end
        )

        if first_line is None:
            return None
        else:
            cpg_sites = self.get_cpg_sites(
                bedgraph_fname, first_line, last_line)
        return cpg_sites

    def bisection_search(self, bedgraph_fname, chrom_first_line, chrom_last_line, region_start, region_end):
        '''
        Search through a bedgraph to find the first and last lines contained within
        a given genomic region.
        '''

        # Get the data from the first and last lines of the region to be searched
        first_line_list = linecache.getline(bedgraph_fname, chrom_first_line).strip().split()
        last_line_list = linecache.getline(bedgraph_fname, chrom_last_line).strip().split()

        region_first_line = 0
        region_last_line = 0

        if region_start <= int(first_line_list[1]):
            # Region starts before first CpG in the chromosome
            region_first_line = chrom_first_line
        else:
            found = False
            top_line = chrom_first_line
            bottom_line = chrom_last_line
            mid_line = int(floor((top_line + bottom_line)/2))

            while (bottom_line != top_line + 1) and not found:
                line = linecache.getline(bedgraph_fname, mid_line).strip().split()
                position = int(line[1])

                if position > region_start:
                    # Take top half
                    bottom_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position < region_start:
                    # Take bottom half
                    top_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position == region_start:
                    # Exact hit
                    region_first_line = mid_line
                    found = True

            if not found:
                region_first_line = bottom_line

        if region_end >= int(last_line_list[1]):
            # Region ends after last CpG in the chromosome
            region_last_line = chrom_last_line
        else:
            found = False
            top_line = region_first_line
            bottom_line = chrom_last_line
            mid_line = int(floor((top_line + bottom_line)/2))

            while (bottom_line != top_line + 1) and not found:
                line = linecache.getline(bedgraph_fname, mid_line).strip().split()
                position = int(line[1])

                if position > region_end:
                    # Take top half
                    bottom_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position < region_end:
                    # Take bottom half
                    top_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position == region_end:
                    region_last_line = mid_line
                    found = True

            if not found:
                region_last_line = top_line

        if region_first_line == region_last_line:
            # Either there are no CpGs in this region, or this line contains the
            # only CpG in the region
            line = linecache.getline(bedgraph_fname, region_first_line).strip().split()
            position = int(line[1])
            if region_start <= position <= region_end:
                # Only CpG in the region
                pass
            else:
                # No CpGs in the region
                return None,None

        return [region_first_line, region_last_line]

    def get_cpg_sites(self, bedgraph_fname, first_line, last_line):
        cpg_sites = []
        for i in range(first_line, last_line+1):
            line_list = linecache.getline(bedgraph_fname, i).strip().split()
            coord = int(line_list[1])
            meth_level = float(line_list[3])
            cpg_sites.append((coord,meth_level))
        return cpg_sites

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

def parse_gtf(gtf_fname):
    logging.info("Parsing GTF file")

    gtf_dict = {}
    with open(gtf_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split("\t")
            if line_list[2] not in ["transcript", "exon"]:
                continue
            try:
                attrs = get_attr_dict(line_list[-1])
            except ValueError:
                print(tb.format_exc())
                print(line_list[-1])
                sys.exit()

            if line_list[2] == "transcript":
                chromosome = line_list[0]
                if "chr" not in chromosome:
                    chromosome = "chr" + chromosome

                gtf_dict[attrs["transcript_id"]] = {
                    "chromosome":chromosome,
                    "start":int(line_list[3]),
                    "end":int(line_list[4]),
                    "strand":line_list[6],
                    "exons":[],
                    "gene_name":attrs["gene_name"]
                }
            elif line_list[2] == "exon":
                gtf_dict[attrs["transcript_id"]]["exons"]\
                    .append((int(line_list[3]), int(line_list[4])))

    logging.info("Finished parsing GTF file")

    return gtf_dict

def get_bedgraph_dict(bedgraph_list, group_list, bedgraph_labels):
    logging.info("Creating bedgraph dictionary")

    bedgraph_dict = {l:[] for l in bedgraph_labels}

    for g,b,l in zip(group_list, bedgraph_list, bedgraph_labels):
        bedgraph_dict[l] = (b,g)

    logging.info("Finished bedgraph dictionary")
    return bedgraph_dict

def parse_idx(idx_fname):
    logging.info("Fetching bedgraph index info from {}".format(idx_fname))

    index_dict = {}
    with open(idx_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split(",")
            index_dict[line_list[0]] = {
                "first":int(line_list[1]),
                "last":int(line_list[2])
            }

    logging.info("Finished")
    return index_dict

def get_bedgraph_idxs(bedgraph_dict):
    bedgraph_idxs = {}
    for label, (bedgraph_fname, group) in bedgraph_dict.iteritems():
        idx_fname = bedgraph_fname + ".idx"
        bedgraph_idxs[label] = parse_idx(idx_fname)
    return bedgraph_idxs

def get_meth_data(transcript_list_fname, gtf_dict, add_left, add_right, bedgraph_dict, num_procs):
    logging.info("Starting to get methylation data")

    logging.info("Reading in transcript list")
    with open(transcript_list_fname, 'r') as f:
        transcripts = [l.strip() for l in f]

    logging.info("Creating transcript dictionary")
    transcript_dict = {t:gtf_dict[t] for t in transcripts}

    logging.info("Creating bedgraph index dictionaries")
    bedgraph_idxs = get_bedgraph_idxs(bedgraph_dict)

    job_queue = mp.Queue()
    results_queue = mp.Queue()

    logging.info("Starting methylation data retrieval processes")

    procs = [GetMethProcess(job_queue, results_queue) for _ in range(num_procs)]
    for p in procs:
        p.start()

    logging.info("All procecesses started")

    logging.info("Starting to put jobs on queue")

    for t in transcripts:
        trans_info = transcript_dict[t].copy()
        transcript_dict[t]["methylation"] = {}
        for l,(bg, group) in bedgraph_dict.iteritems():
            idx = bedgraph_idxs[l]
            job_queue.put([t, l, trans_info, add_left, add_right, bg, idx])
            transcript_dict[t]["methylation"][l] = {"group":group, "cpgs":[]}

    logging.info("All jobs on queue")

    for _ in range(num_procs):
        job_queue.put("STOP")

    logging.info("Retrieving results from queue")

    expected_results = len(transcripts)*len(bedgraph_dict)
    count = 0
    while count < expected_results:
        t, l, cpgs = results_queue.get()
        transcript_dict[t]["methylation"][l]["cpgs"] = cpgs
        count += 1

    logging.info("All results retrieved")

    return transcript_dict

def plot_methylation(transcript, info, add_left, add_right, color_dict, img_fname):
    title = "{} ({})".format(transcript, info["gene_name"])

    fig = plt.figure()
    fig.set_size_inches(12,6)
    plt.suptitle(title)

    gs = gridspec.GridSpec(2, 1, height_ratios = [8, 1])
    line_ax = plt.subplot(gs[0])
    transcript_ax = plt.subplot(gs[1], sharex=line_ax)

    xmin = info["start"] - add_left
    xmax = info["end"] + add_right

    line_ax.set_xlim([xmin, xmax])
    line_ax.set_ylim([0,105])

    line_ax.set_ylabel("Methylation")

    for l,d in info["methylation"].iteritems():
        color = color_dict[d["group"]]
        data = d["cpgs"]
        x = []
        y = []
        for p in data:
            x.append(p[0])
            y.append(p[1])
        line_ax.plot(x, y, '-', color=color, ms=0, lw=1, alpha=0.5)

    line_ax.get_xaxis().set_visible(False)
    transcript_ax.get_yaxis().set_visible(False)

    for s in transcript_ax.spines.values():
        s.set_color('none')

    transcript_ax.set_ylim([0,1])

    transcript_start = info["start"]
    transcript_end = info["end"]

    total_length = (transcript_end + add_right) - (transcript_start - add_left)

    transcript_ax.axhline(y=0.5, color='black', lw=1, zorder=0)
    print(float(add_left)/total_length)
    print(float(transcript_end - transcript_start + add_left)/total_length)
    transcript_ax.axhline(
        y=0.5,
        xmin=float(add_left)/total_length,
        xmax=float(transcript_end - transcript_start + add_left)/total_length,
        color='green',
        lw=3,
        zorder=0
    )

    transcript_ax.set_xlabel("Length (bp)\n{}:{}-{}".format(info["chromosome"], transcript_start, transcript_end))
    transcript_ax.set_xticklabels(np.arange(0, total_length, 1000))

    for exon_start, exon_end in info["exons"]:
        print(exon_start)
        print(exon_end)
        exon_length = exon_end - exon_start
        rect = mpatches.Rectangle((exon_start, 0), exon_length, 1, linewidth=0, color='green')
        transcript_ax.add_patch(rect)

        highlight = mpatches.Rectangle((exon_start, 0), exon_length, 105, linewidth=0, alpha=0.2, color='green')
        line_ax.add_patch(highlight)

    plt.savefig(img_fname)

    plt.close()

def plot_all_transcripts(transcript_dict, add_left, add_right, color_dict, output_dir):
    for transcript, info in transcript_dict.iteritems():
        img_fname = os.path.join(output_dir, transcript+".svg")
        plot_methylation(transcript, info, add_left, add_right, color_dict, img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("transcript_list_fname")
    parser.add_argument("gtf_fname")
    parser.add_argument("num_procs", type=int)
    parser.add_argument("-b", "--bedgraph", action="append")
    parser.add_argument("-g", "--group", action="append")
    parser.add_argument("-L", "--bedgraphLabel", action="append")
    parser.add_argument("-l", "--addLeft", type=int, default=1000)
    parser.add_argument("-r", "--addRight", type=int, default=1000)
    parser.add_argument("output_dir")
    args = parser.parse_args()

    gtf_dict = parse_gtf(args.gtf_fname)

    color_dict = {"B":"blue", "T":"red"}

    bedgraph_dict = get_bedgraph_dict(args.bedgraph, args.group, args.bedgraphLabel)

    transcript_dict = get_meth_data(args.transcript_list_fname, gtf_dict, args.addLeft, args.addRight, bedgraph_dict, args.num_procs)

    plot_all_transcripts(transcript_dict, args.addLeft, args.addRight, color_dict, args.output_dir)
