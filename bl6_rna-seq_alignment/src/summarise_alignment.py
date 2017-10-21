'''
Script to summarise read mapping with STAR to produce two barcharts.
One barchart showing percentages, one showing number of reads.
Each sample will show number unmapped, number assigned to rRNA, number uniquely
mapped, and number mapping to multiple locations.
Also produce density plots of alignment scores
'''

import sys
import os
import re
import json
import argparse
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def parse_star_log(log_fname):
    '''
    Parse a STAR Log.final.out file to obtain total reads and percentage
    of reads unmapped
    '''
    with open(log_fname, 'r') as f:
        for i,line in enumerate(f):
            line = line.strip()
            if i == 5:
                total = int(re.split("\s+", line)[-1])
            elif i == 25:
                too_many_hits = int(re.split("\s+", line)[-1])
            elif i == 28:
                too_many_mismatches_pc = float(re.split("\s+", line)[-1][0:-1])
            elif i == 29:
                too_short_pc = float(re.split("\s+", line)[-1][0:-1])
            elif i == 30:
                other_pc = float(re.split("\s+", line)[-1][0:-1])


    unmapped_pc = ((100.0*too_many_hits)/total) \
                    + too_many_mismatches_pc \
                    + too_short_pc \
                    + other_pc

    return total, unmapped_pc

def parse_num_in(num_in_fname):
    with open(num_in_fname, 'r') as f:
        return int(f.readline().strip())

def parse_hits_summary(hits_summ_fname):
    with open(hits_summ_fname, 'r') as f:
        f.readline()
        multimapped = 0
        uniq = 0
        for line in f:
            hits, num_reads = [int(x) for x in line.strip().split(",")]
            if hits == 1:
                uniq = num_reads
            else:
                multimapped += num_reads
    return multimapped, uniq

def get_data_dict(star_dir, rrna_dir):
    '''
    Produce a dictionary with all of the data required for plotting, structured
    as follows:
    KEY "read_nums"
        -> KEY <sample name>
            -> KEY "total"
                -> VALUE <total reads>
            -> KEY "unmapped"
                -> VALUE <num read unmapped>
            -> KEY  "rRNA"
                -> VALUE <num rRNA reads>
            -> KEY "multimapped"
                -> VALUE <num reads multimapped>
            -> KEY "uniq"
                -> VALUE <num reads uniquely mapped
    KEY "read_percents"
        As above, without total
    KEY "alignment_scores"
        -> KEY <sample name>
            -> KEY <score>
                -> VALUE <count>
    '''
    samples = os.listdir(star_dir)
    data_dict = {
        "read_nums":{
            s:{
                k:0 for k in ["total", "unmapped", "rRNA", "multimapped", "uniq"]
            } for s in samples
        },
        "read_percents":{
            s:{
                k:0 for k in ["unmapped", "rRNA", "multimapped", "uniq"]
            } for s in samples
        },
        "alignment_scores":{s:{} for s in samples}
    }

    # First parse STAR's Log.final.out files to extract the total number of
    # reads and percentage unmapped
    for s in samples:
        log_fname = os.path.join(star_dir, s, "Log.final.out")
        total, unmapped_pc = parse_star_log(log_fname)
        data_dict["read_nums"][s]["total"] = total
        data_dict["read_percents"][s]["unmapped"] = unmapped_pc
        data_dict["read_nums"][s]["unmapped"] = int((unmapped_pc/100.0)*total)

    # Now parse the sumamry files for the RSeQC output
    for s in samples:
        num_in_fname = os.path.join(rrna_dir, s, "number_in.txt")
        hits_summ_fname = os.path.join(rrna_dir, s, "ex_hits_summary.csv")

        rrna = parse_num_in(num_in_fname)
        multimapped, uniq = parse_hits_summary(hits_summ_fname)

        data_dict["read_nums"][s]["rRNA"] = rrna
        data_dict["read_nums"][s]["multimapped"] = multimapped
        data_dict["read_nums"][s]["uniq"] = uniq

        total = data_dict["read_nums"][s]["total"]

        data_dict["read_percents"][s]["rRNA"] = (100.0*rrna)/total
        data_dict["read_percents"][s]["multimapped"] = (100.0*multimapped)/total
        data_dict["read_percents"][s]["uniq"] = (100.0*uniq)/total

    return data_dict

def hbar_percentages(pc_data, base_plot_fname):
    '''
    Draw a horizontal barchart showing the breakdown of read mapping for each
    sample, as a percentage
    '''

    fig, ax = plt.subplots(figsize = (8,10))
    plt.title("RNA-seq Alignment Breakdown (%)")

    for i,sample in enumerate(sorted(pc_data.keys())):
        unmapped = pc_data[sample]["unmapped"]
        rrna = pc_data[sample]["rRNA"]
        multimapped = pc_data[sample]["multimapped"]
        uniq = pc_data[sample]["uniq"]

        left = 0.0
        for d,col,a in zip(
                [uniq, multimapped, rrna, unmapped],
                ["green", "blue", "red", "grey"],
                [0.7,0.4,0.4,0.4]):
            ax.barh(i, d, left = left, color = col, alpha = a, align = 'center')
            left += d

    ax.set_yticks(range(len(pc_data.keys())))
    ax.set_yticklabels(sorted(pc_data.keys()))

    ax.set_xlim([0,100])
    ax.set_ylim([-1, len(pc_data.keys())])

    ax.set_ylabel("Sample")
    ax.set_xlabel("Percent of reads")

    legend_patches = [
                        mpatches.Patch(color = 'green', alpha = 0.7, label = "Uniquely mapped"),
                        mpatches.Patch(color = 'blue', alpha = 0.7, label = "Multimapped"),
                        mpatches.Patch(color = 'red', alpha = 0.7, label = "rRNA"),
                        mpatches.Patch(color = 'grey', alpha = 0.7, label = "Unmapped")
    ]
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    print("Done, saving figure")
    sys.stdout.flush()

    plot_fname = "pc_aln_breakdown_"+base_plot_fname
    plt.savefig(plot_fname, bbox_inches='tight')

def hbar_read_nums(num_data, base_plot_fname):
    '''
    Draw a horizontal barchart showing the breakdown of read mapping for each
    sample, as numbers of raw reads
    '''

    fig, ax = plt.subplots(figsize=(8,10))
    plt.title("RNA-seq Alignment Breakdown")

    for i,sample in enumerate(sorted(num_data.keys())):
        unmapped = num_data[sample]["unmapped"]
        rrna = num_data[sample]["rRNA"]
        multimapped = num_data[sample]["multimapped"]
        uniq = num_data[sample]["uniq"]

        left = 0.0
        for d,col,a in zip(
                [uniq, multimapped, rrna, unmapped],
                ["green", "blue", "red", "grey"],
                [0.7,0.4,0.4,0.4]):
            ax.barh(i, d, left=left, color=col, alpha=a, align='center')
            left += d

    ax.set_yticks(range(len(num_data.keys())))
    ax.set_yticklabels(sorted(num_data.keys()))

    ax.set_ylim([-1, len(num_data.keys())])

    ax.set_ylabel("Sample")
    ax.set_xlabel("Number of reads")

    legend_patches = [
                        mpatches.Patch(color='green', alpha=0.7, label="Uniquely mapped"),
                        mpatches.Patch(color='blue', alpha=0.7, label="Multimapped"),
                        mpatches.Patch(color='red', alpha=0.7, label="rRNA"),
                        mpatches.Patch(color='grey', alpha=0.7, label="Unmapped")
    ]
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    print("Done, saving figure")
    sys.stdout.flush()

    plot_fname = "rawnum_aln_breakdown_"+base_plot_fname
    plt.savefig(plot_fname, bbox_inches = 'tight')
'''
def bar_align_scores(aln_score_data, base_plot_fname):

    Draw a barchart showing distribution of the alignment scores of the non-rRNA
    reads for each sample

    fig, axarr = plt.subplots(4, 3, figsize = (8,8), sharex = True, sharey = True)

    plt.title("Alignment Score Distributions")

    for i,sample in enumerate(sorted(aln_score_data.keys())):
        ax = axarr[int(i/3), int(i%3)]
        inds = [int(x) for x in sorted(aln_score_data[sample].keys())]
        heights = [aln_score_data[sample][k] for k in inds]

        ax.bar(inds, heights, width = 1, color = 'blue', lw = 0)
        #ax.bar(range(len(heights)), heights, width = 1, color = 'blue', lw = 0)

        ax.set_title(sample)

    for i in [0,1,2]:
        axarr[3,i].set_xlabel("Alignment Score")

    for i in [0,1,2,3]:
        axarr[i,0].set_ylabel("Number of Alignments")

    ax.set_xticks

    plot_fname = "aln_scores_" + base_plot_fname
    plt.savefig(plot_fname, bbox_inches = 'tight')
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_json")
    parser.add_argument("star_dir")
    parser.add_argument("rrna_dir")
    parser.add_argument("plot_fname")
    args = parser.parse_args()

    try:
        with open(args.data_json, 'r') as f:
            data_dict = json.load(f)
    except IOError:
        # Need to get the data ourselves
        data_dict = get_data_dict(args.star_dir, args.rrna_dir)
        with open(args.data_json, 'wa') as f:
            json.dump(data_dict, f)

    hbar_percentages(data_dict["read_percents"], args.plot_fname)
    hbar_read_nums(data_dict["read_nums"], args.plot_fname)
    #bar_align_scores(data_dict["alignment_scores"], args.plot_fname)
