'''
Plot summaries of reconstructed transcriptomes, including:
- transcripts per gene
- exons per transcripts
- transcripts on each chromosome
- number of genes, transcripts, exons
Include a reference transcriptome as a comparison. Needs a JSON file as input,
which should contain a dictionary as follows:
{file_label:{"fname":"/path/to/gtf", "color":"rgb_color_code"}}
'''

import re
import json
import subprocess
from StringIO import StringIO
import logging

from numpy import log10, array
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

logging.basicConfig(level=logging.INFO)

def get_gtf_summ(gtf_fname):
    logging.info("Getting GTF summary for file {}".format(gtf_fname))

    gtf_summ = {
        "transcripts_per_gene":{},
        "exons_per_transcript":{},
        "transcripts_per_chromosome":{},
        "total_exons":-1,
        "total_transcripts":-1,
        "total_genes":-1}

    gtf_summ["total_exons"] = int(
        subprocess.check_output(
            """cat %s |
            sed '/^#/d' |
            grep -P '\sexon\s' |
            awk '{print $1,$4,$5,$7}' |
            sort -k1,1 -k2,2n | uniq | wc -l""" % gtf_fname,
            shell=True))

    gtf_summ["total_transcripts"] = int(
        subprocess.check_output(
            """cat %s |
            sed '/^#/d' |
            grep -P '\stranscript\s' |
            awk '{print $12}' |
            sort | uniq | wc -l""" % gtf_fname,
            shell=True))

    gtf_summ["total_genes"] = int(
        subprocess.check_output(
            """cat %s |
            sed '/^#/d' |
            grep -P '\stranscript\s' |
            awk '{print $10}' |
            sort | uniq | wc -l""" % gtf_fname,
            shell=True))

    gtf_summ["transcripts_per_chromosome"] = {
        re.sub("chr", "", chrom):int(count)
        for chrom, count in map(
            lambda l: l.strip().split(),
            StringIO(
                subprocess.check_output(
                    """cat %s |
                    sed '/^#/d' |
                    grep -P '\stranscript\s' |
                    awk '{print $1,$12}' |
                    sort | uniq |
                    awk '{print $1}' |
                    sort | uniq -c |
                    sed '/[JMGS]/d' |
                    awk '{print $2,$1}'""" % gtf_fname,
                    shell=True)))}

    gtf_summ["exons_per_transcript"] = {
        int(number):int(frequency)
        for number,frequency in map(
            lambda l: l.strip().split(),
            StringIO(
                subprocess.check_output(
                    """cat %s |
                    sed '/^#/d' |
                    grep -P '\sexon\s' |
                    awk '{print $12}' |
                    sort | uniq -c |
                    awk '{print $1}' |
                    sort | uniq -c |
                    awk '{print $2,$1}'""" % gtf_fname,
                    shell=True)))}

    gtf_summ["transcripts_per_gene"] = {
        int(number):int(frequency)
        for number,frequency in map(
            lambda l: l.strip().split(),
            StringIO(
                subprocess.check_output(
                    """cat %s |
                    sed '/^#/d' |
                    grep -P '\stranscript\s' |
                    awk '{print $10}' |
                    sort | uniq -c |
                    awk '{print $1}' |
                    sort | uniq -c |
                    awk '{print $2,$1}'""" % gtf_fname,
                    shell=True)))}

    return gtf_summ

def get_plot_data(input_dict):
    raw_data = {label:get_gtf_summ(d["fname"]) for label,d in input_dict.items()}

    #TODO: Make it work for Ensembl GTFs?
    #raw_data[ref_label] = {"color":"black", "data":get_gtf_summ(ref_gtf)}

    plot_data = {
        "transcripts_per_gene":{},
        "exons_per_transcript":{},
        "transcripts_per_chromosome":{},
        "total_exons":{},
        "total_transcripts":{},
        "total_genes":{}}

    for label,d in raw_data.iteritems():
        plot_data["transcripts_per_gene"][label] = d["transcripts_per_gene"]
        plot_data["exons_per_transcript"][label] = d["exons_per_transcript"]
        plot_data["transcripts_per_chromosome"][label] = d["transcripts_per_chromosome"]
        for plot_type in ["total_genes", "total_transcripts", "total_exons"]:
            plot_data[plot_type].update({label:d[plot_type]})

    return plot_data

def cmp_chrs(c1, c2):
    try:
        c1 = int(c1)
    except ValueError:
        c1 = ord(c1)

    try:
        c2 = int(c2)
    except ValueError:
        c2 = ord(c2)

    return c1-c2

def draw_plots(plot_data, input_dict, log, img_fname):
    font = {'size':6}
    matplotlib.rc('font', **font)

    fig = plt.figure()
    fig.set_size_inches(12,6)
    gs = gridspec.GridSpec(2, 3)
    gs.update(wspace=0.3, hspace=0.3)

    ax_dict = {
        "transcripts_per_gene":plt.subplot(gs[0,0]),
        "exons_per_transcript":plt.subplot(gs[0,1]),
        "transcripts_per_chromosome":plt.subplot(gs[0,2]),
        "total_exons":plt.subplot(gs[1,2]),
        "total_transcripts":plt.subplot(gs[1,1]),
        "total_genes":plt.subplot(gs[1,0])}

    title_dict = {
        "transcripts_per_gene":"Transcripts per gene",
        "exons_per_transcript":"Exons per transcript",
        "transcripts_per_chromosome":"Transcripts per chromosome",
        "total_exons":"Total exons",
        "total_transcripts":"Total transcripts",
        "total_genes":"Total genes"}

    xlabel_dict = {
        "transcripts_per_gene":"# transcripts per gene",
        "exons_per_transcript":"# exons per transcript",
        "transcripts_per_chromosome":"Chromosome",
        "total_exons":"Sample",
        "total_transcripts":"Sample",
        "total_genes":"Sample"}

    if log:
        ylabel_dict = {
            "transcripts_per_gene":"$\log_{10}(frequency)$",
            "exons_per_transcript":"$\log_{10}(frequency)$",
            "transcripts_per_chromosome":"$\log_{10}($Number of transcripts$)$",
            "total_exons":"# exons",
            "total_transcripts":"# transcripts",
            "total_genes":"# genes"}
    else:
        ylabel_dict = {
            "transcripts_per_gene":"Frequency",
            "exons_per_transcript":"Frequency",
            "transcripts_per_chromosome":"Frequency",
            "total_exons":"# exons",
            "total_transcripts":"# transcripts",
            "total_genes":"# genes"}

    for plot_type in ["transcripts_per_gene", "exons_per_transcript"]:
        for label, d in plot_data[plot_type].iteritems():
            x = sorted(d.keys())
            if log:
                y = [log10(d[i]) for i in x]
            else:
                y = [d[i] for i in x]
            ax_dict[plot_type].plot(
                x, y, '-o',
                color=input_dict[label]['color'],
                alpha=0.7,
                markersize=2,
                lw=1)
            ax_dict[plot_type].set_title(title_dict[plot_type])
            ax_dict[plot_type].set_xlabel(xlabel_dict[plot_type])
            ax_dict[plot_type].set_ylabel(ylabel_dict[plot_type])

    for label, d in plot_data["transcripts_per_chromosome"].iteritems():
        chroms = sorted(d.keys(), cmp=cmp_chrs)
        x = range(len(chroms))
        if log:
            y = [log10(d[c]) for c in chroms]
        else:
            y = [d[c] for c in chroms]
        ax_dict["transcripts_per_chromosome"].plot(
            x, y, '-o',
            color=input_dict[label]['color'],
            alpha=0.7,
            markersize=2,
            lw=1,
            label=label)
        ax_dict["transcripts_per_chromosome"].set_xticks(x)
        ax_dict["transcripts_per_chromosome"].set_xticklabels(chroms)
        ax_dict["transcripts_per_chromosome"].set_title(title_dict["transcripts_per_chromosome"])
        ax_dict["transcripts_per_chromosome"].set_xlabel(xlabel_dict["transcripts_per_chromosome"])
        ax_dict["transcripts_per_chromosome"].set_ylabel(ylabel_dict["transcripts_per_chromosome"])

    width = 1.0
    for plot_type in ["total_genes", "total_transcripts", "total_exons"]:
        data = plot_data[plot_type]
        labels = sorted(data.keys())
        idxs = array(range(len(labels)))
        heights = [data[l] for l in labels]
        bars = ax_dict[plot_type].bar(
            idxs+width, heights, width, color='blue', alpha=0.7, lw=0)
        for i,b in enumerate(bars):
            b.set_color(input_dict[labels[i]]['color'])
        ax_dict[plot_type].set_xticks(idxs + width)
        ax_dict[plot_type].set_xticklabels(labels, rotation=60, ha='right')
        ax_dict[plot_type].set_title(title_dict[plot_type])
        ax_dict[plot_type].set_xlabel(xlabel_dict[plot_type])
        ax_dict[plot_type].set_ylabel(ylabel_dict[plot_type])

    legend_patches = [
        mpatches.Patch(
            label=l,
            color=input_dict[l]['color']
        ) for l in sorted(input_dict.keys())]

    ax_dict["transcripts_per_chromosome"].legend(
        handles=legend_patches,
        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=8)

    plt.savefig(img_fname, bbox_inches='tight')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json")
    #parser.add_argument("ref_gtf")
    #parser.add_argument("-r", "--ref_label", type=str, default="REF")
    parser.add_argument("--log", action="store_true")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    with open(args.input_json, 'r') as f:
        input_dict = json.load(f)

    plot_data = get_plot_data(input_dict)

    draw_plots(plot_data, input_dict, args.log, args.img_fname)
