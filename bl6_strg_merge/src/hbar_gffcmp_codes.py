'''
Given a set of labelled gffcompare tmap files, draw a horizontal barchart
showing the breakdown of classification codes for each.
'''

import json
import logging
import subprocess
from StringIO import StringIO

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_prop_dict(tmap_fname):
    logging.info("Getting data from file {}".format(tmap_fname))
    count_dict = {
        code:int(count)
        for code,count in map(
            lambda l: l.strip().split(),
            StringIO(
                subprocess.check_output(
                    """cat %s |
                    sed '1d' |
                    awk '{print $3}' |
                    sort | uniq -c |
                    awk '{print $2,$1}'""" % tmap_fname,
                    shell=True)))}

    total = sum(count_dict.values())
    prop_dict = {c:100.0*float(x)/total for c,x in count_dict.items()}
    return prop_dict, total

def get_plot_data(input_dict, classes):
    plot_data = {}
    totals = {}

    for l,fn in input_dict.items():
        data,tot = get_prop_dict(fn)
        plot_data[l] = data
        totals[l] = tot

    return plot_data, totals

def draw_hbars(plot_data, totals, label_list, classes, img_fname):
    matplotlib.rc("font", **{"size":16})
    cmap = matplotlib.cm.get_cmap('Greens')
    color_dict = {}
    for i,c in enumerate(["e", "j", "c", "="]):
        color_dict[c] = cmap(float(i+1)/5)

    color_dict.update({
        "i":"blue",
        "x":"purple",
        "u":"red",
        "o":"orange",
        "p":"gray",
        "s":"beige"
    })

    explanation_dict = {
        "=":"complete match",
        "c":"query contained in ref",
        "j":"possible novel isoform",
        "e":"possible pre-mRNA",
        "i":"intronic",
        "o":"exon overlap",
        "p":"possible polymerase run-on",
        "u":"novel transcript",
        "x":"antisense exon overlap",
        "s":"probable false positive"
    }

    num_bars = len(label_list)

    fig = plt.figure(figsize = [12,1*num_bars])

    axis = fig.add_axes([0.05,0.05,0.9,0.9])

    axis.set_ylabel("Region set")
    axis.set_xlabel("Percentage of transcripts")

    height = 0.8

    axis.set_xlim([-5,105])
    axis.set_ylim([0.5, num_bars+1])

    axis.set_yticks([i + 1 for i in range(num_bars)])
    ytick_labels = [
        "{} ({})".format(label, totals[label]) for label in reversed(label_list)]
    axis.set_yticklabels(ytick_labels)

    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    for i,label in enumerate(reversed(label_list)):
        data = plot_data[label]
        sum_width = 0
        for c in classes:
            percentage = data[c]
            axis.barh(
                i+1, percentage,
                height=height, left=sum_width,
                color=color_dict[c], alpha=0.7, lw=0)

            if percentage > 5:
                axis.text(
                    sum_width + percentage/2, i + 1,
                    "%s\n%.2f" % (c, percentage),
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=12)

            sum_width += percentage


    legend_patches = [
        mpatches.Patch(
            color=color_dict[c],
            label="{}: {}".format(c, explanation_dict[c]))
        for c in classes]

    plt.legend(
        handles=legend_patches,
        loc=2,
        bbox_to_anchor=(1.05, 1.05),
        prop={'size':16})

    plt.savefig(img_fname, bbox_inches='tight')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    classes = ["=", "c", "j", "e", "i", "x", "u", "o", "p", "s"]

    with open(args.input_json, 'r') as f:
        input_dict = json.load(f)

    plot_data, totals = get_plot_data(input_dict, classes)

    draw_hbars(plot_data, totals, args.label, classes, args.img_fname)
