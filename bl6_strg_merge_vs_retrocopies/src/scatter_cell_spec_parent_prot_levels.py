'''
Input:
- sets of B-specific, T-specific, and cell-shared transcripts
- proteome data with gene names
- mapping between transcript IDs and gene names
Output:
- scatter plot of protein levels for each set of transcripts

Assume we are getting median-normalised, log2 transformed proteome data
'''

import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# tested: yes
def get_attr_dict(attr_str):
    return dict(
        map(
            lambda y: re.split(' *"', y)[0:-1],
            map(
                lambda x: x.strip(),
                attr_str.split(";")[0:-1])))

# tested: yes
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

# tested: yes
def get_protein_data(proteome_tsv_fname):
    protein_data = {}
    with open(proteome_tsv_fname, 'r') as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            line_list = line.strip().split("\t")
            line_dict = {h:v for h,v in zip(headers, line_list)}

            #if line_dict["GN"] == "NULL" \
            #or line_dict["OS"] != "Mus_musculus" \
            #or "(Fragment)" in line_dict["Description"]:
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

# tested: yes
def get_transcript_data(transcript_list, transcript2gene, protein_data):
    transcript_data = {
        sex:{cell:[] for cell in ["B", "T"]}
        for sex in ["male", "female"]}

    with open(transcript_list, 'r') as f:
        for line in f:
            transcript = line.strip()
            try:
                gene_name = transcript2gene[transcript]
            except KeyError:
                continue

            try:
                prot_vals = protein_data[gene_name]
            except KeyError:
                prot_vals = {
                    sex:{cell:[0.0] for cell in ["B", "T"]}
                    for sex in ["male", "female"]}

            for sex in ["male", "female"]:
                for cell in ["B", "T"]:
                    transcript_data[sex][cell] += prot_vals[sex][cell]

    return transcript_data

# tested: yes
def get_plot_data(B_spec, T_spec, cell_shared, transcript2gene, protein_data):
    plot_data = {
        "B_spec":get_transcript_data(B_spec, transcript2gene, protein_data),
        "T_spec":get_transcript_data(T_spec, transcript2gene, protein_data),
        "cell_shared":get_transcript_data(cell_shared, transcript2gene, protein_data)
    }

    return plot_data

# tested: no
def draw_scatter(plot_data, img_fname):
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(10,6)

    color_dict = {"B_spec":"blue", "T_spec":"red", "cell_shared":"gray"}
    size_dict = {"B_spec":8, "T_spec":8, "cell_shared":2}
    marker_dict = {"male":"^", "female":"o"}

    for label, transcript_data in plot_data.iteritems():
        for sex, d in transcript_data.iteritems():
            x = d["B"]
            y = d["T"]
            ax.scatter(
                x, y,
                color=color_dict[label], alpha=0.6,
                marker=marker_dict[sex], s=size_dict[label],
                lw=0, label=label+" "+sex)

    lim = max(ax.get_xlim()+ax.get_ylim())
    ax.set_xlim([-0.5,lim])
    ax.set_ylim([-0.5,lim])
    ax.plot([0,lim], [0,lim], color='k', linestyle='dashed', lw=0.5)

    ax.set_xlabel("B cell $\log_2(abundance)$")
    ax.set_ylabel("T cell $\log_2(abundance)$")
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("B_spec")
    parser.add_argument("T_spec")
    parser.add_argument("cell_shared")
    parser.add_argument("gtf_fname")
    parser.add_argument("proteome_tsv_fname")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    transcript2gene = get_transcript_to_gene_dict(args.gtf_fname)

    protein_data = get_protein_data(args.proteome_tsv_fname)

    plot_data = get_plot_data(
        args.B_spec, args.T_spec, args.cell_shared,
        transcript2gene, protein_data)

    draw_scatter(plot_data, args.img_fname)
