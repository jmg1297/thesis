'''
Given a TSV file containing protein abundances, sets of transcript IDs, and a
GTF file linking transcript IDs to gene names, plot protein abundance density
for each set of transcripts.
'''

import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import gaussian_kde

def get_attr_dict(attr_str):
    return dict(
        map(
            lambda y: re.split(' *"', y)[0:-1],
            map(
                lambda x: x.strip(),
                attr_str.split(";")[0:-1])))

def get_transcript_to_gene_dict(gtf_fname):
    transcript2gene = {}
    with open(gtf_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split("\t")
            if line_list[2] != "transcript":
                continue

            attr_dict = get_attr_dict(line_list[-1])

            transcript2gene[attr_dict["transcript_id"]] = attr_dict["gene_name"]

    return transcript2gene

def get_protein_data(proteome_tsv_fname):
    protein_data = {}
    with open(proteome_tsv_fname, 'r') as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            line_list = line.strip().split("\t")
            line_dict = {h:v for h,v in zip(headers, line_list)}

            if line_dict["GN"] == "NULL":
                continue
            else:
                if line_dict["GN"] not in protein_data:
                    protein_data[line_dict["GN"]] = {
                        sex:{cell:[] for cell in ["B", "T"]}
                        for sex in ["male", "female"]}

            for key in ["Tmale", "Bmale", "Tfemale", "Bfemale"]:
                if "B" in key:
                    cell = "B"
                else:
                    cell = "T"

                if "female" in key:
                    sex = "female"
                else:
                    sex = "male"

                try:
                    value = float(line_dict[key])
                except ValueError:
                    value = 0.0

                protein_data[line_dict["GN"]][sex][cell].append(value)

    return protein_data

def get_plot_data(transcript_lists, labels, transcript2gene, protein_data):
    fname_dict = {l:fn for l,fn in zip(labels, transcript_lists)}

    plot_data = {
        s:{
            c:{
                l:[] for l in labels}
            for c in ["B", "T"]}
        for s in ["male", "female"]
    }

    for l,fn in fname_dict.items():
        with open(fn, 'r') as f:
            for line in f:
                transcript = line.strip()
                gene_name = transcript2gene[transcript]
                try:
                    pdata = protein_data[gene_name]
                except KeyError:
                    continue
                for sex, d in pdata.items():
                    for cell, values in d.items():
                        plot_data[sex][cell][l] += values

    for s in ["male", "female"]:
        for c in ["B", "T"]:
            plot_data[s][c]["ALL"] = reduce(
                lambda a,x: a+x,
                [protein_data[k][s][c] for k in protein_data.keys()],
                [])

    return plot_data

def flatten(nested_dict):
    '''Given a nested dictionary where the lowest-level values are lists,
    concatenate all of these lists and return it'''

    if type(nested_dict) == list:
        return nested_dict
    else:
        return reduce(lambda a,y: a+y, map(flatten, nested_dict.values()), [])

def draw_density_plots(plot_data, label_list, colors, img_fname, bandwidth):

    color_dict = {l:c for l,c in zip(label_list, colors)}
    color_dict["ALL"] = 'black'
    label_list.append("ALL")

    fig, axarr = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.set_size_inches(14,8)

    ax_dict = {
        "female":{"B":axarr[0,0], "T":axarr[0,1]},
        "male":{"B":axarr[1,0], "T":axarr[1,1]}
    }

    all_vals = flatten(plot_data)
    xvals = np.linspace(min(all_vals), max(all_vals), 1000)

    for sex, d in plot_data.iteritems():
        for cell, e in d.iteritems():
            for l, data in e.iteritems():
                density = gaussian_kde(data)
                if bandwidth > 0:
                    density.covariance_factor = lambda : bandwidth
                    density._compute_covariance()
                yvals = density(xvals)
                ax_dict[sex][cell].plot(xvals, yvals, color=color_dict[l], lw=1)
                ax_dict[sex][cell].set_title(sex + " " + cell)

    axarr[0,0].set_ylabel("Density")
    axarr[1,0].set_xlabel("$\log_2($protein abundance$)$")
    axarr[1,0].set_ylabel("Density")
    axarr[1,1].set_xlabel("$\log_2($protein abundance$)$")

    patches = [mpatches.Patch(label=l, color=color_dict[l]) for l in label_list]

    axarr[0,1].legend(
        handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.savefig(img_fname)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("proteome_tsv_fname")
    parser.add_argument("gtf_fname")
    parser.add_argument("img_fname")
    parser.add_argument("-t", "--transcript_list", action="append")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("-c", "--color", action="append")
    parser.add_argument("-b", "--bandwidth", type=float, default=-1.0)
    args = parser.parse_args()

    transcript2gene = get_transcript_to_gene_dict(args.gtf_fname)

    protein_data = get_protein_data(args.proteome_tsv_fname)

    plot_data = get_plot_data(args.transcript_list, args.label, transcript2gene, protein_data)

    draw_density_plots(plot_data, args.label, args.color, args.img_fname, args.bandwidth)
