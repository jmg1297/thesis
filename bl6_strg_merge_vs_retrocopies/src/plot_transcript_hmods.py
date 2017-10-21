'''
For a given set of transcripts, for each transcript in the set:
- find overlapping histone mark peaks from several sets of peaks
- plot the transcript structure and overlapping histone mark peaks
- save all plots in a given directory
'''

import sys
import os
import traceback as tb
import re
import linecache
import logging
import subprocess
import multiprocessing as mp
from Queue import Empty
from math import floor
import tempfile as tf

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

logging.basicConfig(level=logging.INFO)

class IntersectionProcess(mp.Process):
    def __init__(self, job_queue, results_queue, add_left, add_right):
        mp.Process.__init__(self)
        self.job_queue = job_queue
        self.results_queue = results_queue
        self.add_left = add_left
        self.add_right = add_right
        logging.info("Initialised {}".format(self.name))

    def run(self):
        logging.info("{} starting to retrieve jobs from queue".format(self.name))

        while True:
            job = self.job_queue.get()
            if job == "STOP":
                logging.info("Received STOP signal, stopping {}".format(self.name))
                self.results_queue.put("STOP")
                return 0
            else:
                transcript, info, mod, group, fname = job
                logging.info("Retrieved {}/{}/{} job, processing".format(transcript, mod, group))
                overlaps = self.get_overlaps(info, fname)
                logging.info("Finished {}/{}/{} job, putting results on queue".format(transcript, mod, group))
                self.results_queue.put([transcript, mod, group, overlaps])

    def get_overlaps(self, info, fname):
        cmd_str = "bedtools intersect -wao -a %s -b %s | awk '{print $4,$5,$6}'"
        overlaps = []
        with tf.NamedTemporaryFile() as tmp:
            bed_str = "\t".join(map(str, [info["chromosome"], info["start"]-self.add_left, info["end"]+self.add_right])) + "\n"
            tmp.write(bed_str)
            tmp.flush()
            raw_overlaps = subprocess.check_output(cmd_str % (tmp.name, fname), shell=True).strip().split("\n")
            for o in raw_overlaps:
                if "." in o:
                    continue
                else:
                    start, end = map(int, o.split()[1:])
                    overlaps.append([start, end])

        return overlaps

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

def get_peak_dict(peak_bed_list, group_list, modification_list):
    # Note the assumption that there is one file per group per modification type
    peak_dict = {
        m:{
            g:"" for g in set(group_list)
        } for m in set(modification_list)
    }

    for b, g, m in zip(peak_bed_list, group_list, modification_list):
        peak_dict[m][g] = b

    return peak_dict

def get_plot_data(transcript_list_fname, gtf_dict, peak_dict, add_left, add_right, num_procs):
    logging.info("Starting to get histone mark data")

    logging.info("Reading in transcript list")
    with open(transcript_list_fname, 'r') as f:
        transcripts = [l.strip() for l in f]

    logging.info("Creating transcript dictionary")
    transcript_dict = {t:gtf_dict[t] for t in transcripts}

    job_queue = mp.Queue()
    results_queue = mp.Queue()

    procs = [IntersectionProcess(job_queue, results_queue, add_left, add_right) for _ in range(num_procs)]
    for p in procs:
        p.start()

    for t in transcripts:
        trans_info = transcript_dict[t].copy()
        transcript_dict[t]["hmods"] = {}
        for m, d in peak_dict.iteritems():
            transcript_dict[t]["hmods"][m] = {}
            for g, fn in d.iteritems():
                job_queue.put([t, trans_info, m, g, fn])
                transcript_dict[t]["hmods"][m][g] = []

    for _ in range(num_procs):
        job_queue.put("STOP")

    stop_count = 0
    while stop_count < num_procs:
        res = results_queue.get()
        if res == "STOP":
            stop_count += 1
        else:
            t, m, g, overlaps = res
            transcript_dict[t]["hmods"][m][g] = overlaps

    return transcript_dict

def plot_overlaps(transcript, info, hmods, groups, add_left, add_right, color_dict, img_fname):
    title = "{} ({})".format(transcript, info["gene_name"])

    total_rows = len(hmods)*len(groups) + 1
    total_length = (info["end"] + add_right) - (info["start"] - add_left)

    fig = plt.figure()
    fig.set_size_inches(12, total_rows + 1)
    plt.suptitle(title)

    gs = gridspec.GridSpec(len(hmods)+1, 1, height_ratios = [2]*len(hmods)+[1])

    transcript_ax = plt.subplot(gs[-1])
    transcript_ax.set_xlim([info["start"] - add_left, info["end"] + add_right])
    transcript_ax.set_xticklabels(np.arange(0, total_length, 1000))
    transcript_ax.set_xlabel("Length (bp)\n{}:{}-{}".format(info["chromosome"], info["start"], info["end"]))
    transcript_ax.get_yaxis().set_visible(False)
    for s in transcript_ax.spines.values():
        s.set_color('none')

    ax_dict = {m:plt.subplot(gs[i], sharex=transcript_ax) for i,m in enumerate(hmods)}

    height_dict = {g:i+0.1 for i,g in enumerate(groups)}

    for m,d in info["hmods"].iteritems():
        ax = ax_dict[m]
        ax.get_xaxis().set_visible(False)
        ax.set_ylim([0,len(groups)])
        ax.set_yticks([i+0.5 for i in range(len(groups))])
        ax.set_yticklabels(groups)
        ax.set_ylabel(m)

        for g,l in d.iteritems():
            height = height_dict[g]
            color = color_dict[g]

            for start, end in l:
                length = end - start
                rect = mpatches.Rectangle((start, height), length, 0.8, linewidth=0, color=color, alpha=0.7)
                ax.add_patch(rect)

    for exon_start, exon_end in info["exons"]:
        print(exon_start)
        print(exon_end)
        exon_length = exon_end - exon_start
        rect = mpatches.Rectangle((exon_start, 0), exon_length, 1, linewidth=0, color='green')
        transcript_ax.add_patch(rect)

        #highlight = mpatches.Rectangle((exon_start, 0), exon_length, 1, linewidth=0, alpha=0.2, color='green')
        for m in hmods:
            highlight = mpatches.Rectangle((exon_start, 0), exon_length, len(groups), linewidth=0, alpha=0.2, color='green')
            ax_dict[m].add_patch(highlight)

    transcript_ax.axhline(y=0.5, color='black', lw=1, zorder=0)
    transcript_ax.axhline(
        y=0.5,
        xmin=float(add_left)/total_length,
        xmax=float(info["end"] - info["start"] + add_left)/total_length,
        color='green',
        lw=3,
        zorder=0
    )
    #arrow = mpatches.Arrow(info["start"], 0.5, 100, 0, width=1000, color='black')
    #transcript_ax.add_patch(arrow)
    if info["strand"] == "+":
        marker = ">"
    elif info["strand"] == "-":
        marker = "<"

    transcript_ax.plot(info["start"], 0.5, marker=marker, ms=10, color="black")
    transcript_ax.plot(info["end"], 0.5, marker=marker, ms=10, color="black")

    plt.savefig(img_fname)

    plt.close()

def plot_all_transcripts(plot_data, color_dict, hmods, groups, add_left, add_right, output_dir):
    for transcript, info in plot_data.iteritems():
        img_fname = os.path.join(output_dir, transcript+".svg")
        plot_overlaps(transcript, info, hmods, groups, add_left, add_right, color_dict, img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("transcript_list_fname")
    parser.add_argument("gtf_fname")
    parser.add_argument("num_procs", type=int)
    parser.add_argument("-b", "--peakBed", action="append")
    parser.add_argument("-g", "--group", action="append")
    parser.add_argument("-m", "--modification", action="append")
    parser.add_argument("-l", "--addLeft", type=int, default=1000)
    parser.add_argument("-r", "--addRight", type=int, default=1000)
    parser.add_argument("output_dir")
    args = parser.parse_args()

    gtf_dict = parse_gtf(args.gtf_fname)

    color_dict = {"B":"blue", "T":"red"}
    hmods = sorted(list(set(args.modification)))
    groups = sorted(list(set(args.group)))

    peak_dict = get_peak_dict(args.peakBed, args.group, args.modification)

    plot_data = get_plot_data(
        args.transcript_list_fname,
        gtf_dict,
        peak_dict,
        args.addLeft,
        args.addRight,
        args.num_procs
    )

    plot_all_transcripts(
        plot_data, color_dict, hmods, groups,
        args.addLeft, args.addRight,
        args.output_dir
    )
