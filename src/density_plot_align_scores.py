'''
Visualise the alignment scores of different sets of retrocopies with their
respective parent transcripts. Alignment score is calculated as % of the
retrocopy that matches the parent, i.e.,
(similarity*match_length)/retrocopy_length
Also show the alignment scores for retrocopies aligned to random transcripts.
'''

import os
import sqlite3
import linecache
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
import numpy as np

def parse_alignment_file(aln_path):
    length_line = linecache.getline(aln_path, 21)
    sim_line = linecache.getline(aln_path, 23)

    length = int(length_line.strip().split()[2])
    similarity = float(re.sub("[\(\)\%\s]", "", sim_line.strip().split()[3]))/100

    return length, similarity


def get_plot_data(
        input_fname_list, label_list,
        alignment_dir, random_alignment_dir, retrocopy_parent_db
    ):
    assert "RANDOM_ALN" not in label_list, "RANDOM_ALN cannot be included as a label"

    plot_data = {l:[] for l in label_list + ["RANDOM_ALN"]}

    aligned_retrocopies = os.listdir(alignment_dir)

    for l, fn in zip(label_list, input_fname_list):
        with open(fn, 'r') as f:
            for line in f:
                line_list = line.strip().split()
                retrocopy_fname = line_list[3] + ".matcher"
                if retrocopy_fname not in aligned_retrocopies:
                    continue
                else:
                    aln_path = os.path.join(alignment_dir, retrocopy_fname)

                length = int(line_list[2]) - int(line_list[1])
                aln_length, similarity = parse_alignment_file(aln_path)

                aln_score = (similarity*aln_length)/length

                plot_data[l].append(aln_score)

    random_alignments = os.listdir(random_alignment_dir)
    random_retrocopies = [".".join(f.split(".")[0:-1]) for f in random_alignments]
    random_aln_paths = [os.path.join(random_alignment_dir, f) for f in random_alignments]

    conn = sqlite3.connect(retrocopy_parent_db)
    cur = conn.cursor()
    rows = cur.execute(
        "SELECT retrogene_id, end_coord-start_coord \
        FROM retrogenes \
        WHERE retrogene_id IN ({})"\
            .format(",".join(["?"]*len(random_retrocopies))),
        random_retrocopies
    )
    length_dict = dict(rows)

    for rand_aln_path, rand_rc in zip(random_aln_paths, random_retrocopies):
        aln_length, similarity = parse_alignment_file(rand_aln_path)
        aln_score = (similarity*aln_length)/length_dict[rand_rc]
        plot_data["RANDOM_ALN"].append(aln_score)

    return plot_data


def draw_density_plot(plot_data, label_list, color_list, plot_title, img_fname):
    assert "black" not in color_list, "Black reserved for random alignments"

    num_plots = len(plot_data)+1
    color_dict = {l:c for l,c in zip(label_list, color_list)}
    color_dict["RANDOM_ALN"] = "black"

    fig, axarr = plt.subplots(num_plots, 1, sharex=True)
    fig.set_size_inches(6,4*(num_plots))

    xvals = np.linspace(0, 1, 1000)

    for label in label_list + ["RANDOM_ALN"]:
        data = plot_data[label]
        density = gaussian_kde(data)
        density.covariance_factor = lambda : 0.2
        density._compute_covariance()
        yvals = density(xvals)
        m = max(yvals)
        yvals = [y/m for y in yvals]
        axarr[0].plot(xvals, yvals, color=color_dict[label], label=label)

    axarr[0].set_ylabel("Normalised density")
    axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    axarr[0].set_ylim([0,1.05])

    for label, ax in zip(label_list + ["RANDOM_ALN"], axarr[1:]):
        data = plot_data[label]
        ax.hist(data, color=color_dict[label], lw=0, bins=50)
        ax.text(x=0.01, y=0.95*ax.get_ylim()[1], s=5, text=label, ha='left', va='top')
        ax.set_ylabel("Frequency")

    axarr[-1].set_xlabel("Alignment score")

    plt.suptitle(plot_title, y=0.92)

    plt.savefig(img_fname, bbox_inches='tight')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFname", action="append", help="BED of retrocopies")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("-c", "--color", action="append")
    parser.add_argument("alignment_dir")
    parser.add_argument("random_alignment_dir")
    parser.add_argument("retrocopy_parent_db")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(
        args.inputFname, args.label,
        args.alignment_dir, args.random_alignment_dir, args.retrocopy_parent_db
    )

    draw_density_plot(
        plot_data, args.label, args.color, args.plot_title, args.img_fname
    )
