'''
Given B and T expressed retrocopies, retrieve the B and T expression values for
each parent transcript. Then draw two scatter plots, one for retrocopy
transcripts sense wrt parent and one antisense wrt parent.
'''

import sqlite3
import json
import subprocess
from collections import OrderedDict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from numpy import log10

def get_plot_data(
        B_spec_mbed, T_spec_mbed, parent_expr_db, strand_combo_dict,
        intrachrom_only, log_tpm
    ):

    rel_directions = ["sense", "antisense"]
    categories = ["B_spec", "T_spec"]
    filenames = [B_spec_mbed, T_spec_mbed]

    plot_data = {rd:{c:OrderedDict() for c in categories} for rd in rel_directions}

    conn = sqlite3.connect(parent_expr_db)
    cur = conn.cursor()

    get_parents_cmd = "cat %s | awk '{print $17}' | sort | uniq | sed -r '/\./d'"
    get_expr_query = "SELECT transcript_id, b_cell, t_cell \
                        FROM transcript_averages \
                        WHERE transcript_id IN ({})"

    for cat, fn in zip(categories, filenames):
        parents = subprocess.check_output(
            get_parents_cmd % fn, shell=True).splitlines()
        rows = cur.execute(get_expr_query.format(",".join(["?"]*len(parents))), parents).fetchall()
        if log_tpm:
            expr_dict = {x[0]:map(lambda z: log10(z+1), x[1:]) for x in rows}
        else:
            expr_dict = {x[0]:x[1:] for x in rows}

        with open(fn, 'r') as f:
            for line in f:
                line_list = line.strip().split()
                if intrachrom_only:
                    if line_list[6] != line_list[13]:
                        continue

                transcript_strand = line_list[5]
                retrocopy_strand = line_list[11]
                parent_strand = line_list[18]
                if parent_strand == ".":
                    continue

                rel_dir = strand_combo_dict[parent_strand][retrocopy_strand][transcript_strand]

                parent = line_list[16]

                plot_data[rel_dir][cat][parent] = expr_dict[parent]

    return plot_data

def find_outliers(X, Y, labels):
    upper_c = 0.01
    lower_c = -upper_c
    step = 0.01
    threshold = 0.9
    total = len(labels)
    #outliers = []

    while True:
        count = 0
        outliers = []
        for x,y,l in zip(X,Y,labels):
            if y > x + lower_c and y < x + upper_c:
                count += 1
            else:
                outliers.append((l, x, y))
        if float(count)/total > threshold:
            print("Finished")
            break
        else:
            upper_c += step
            lower_c = -upper_c

    return upper_c, lower_c, list(set(outliers))

def draw_scatter(plot_data, log_tpm, report_outliers, plot_title, img_fname):
    rel_directions = ["sense", "antisense"]
    categories = ["B_spec", "T_spec"]

    fig, axarr = plt.subplots(2, 1)

    ax_dict = {"sense":axarr[0], "antisense":axarr[1]}

    plt.suptitle(plot_title)

    color_dict = {"B_spec":"blue", "T_spec":"red"}
    label_dict = {
        "B_spec":"B specific \nexpressed \nretrocopy",
        "T_spec":"T specific \nexpressed \nretrocopy"
    }

    for rel_dir,d in plot_data.iteritems():
        all_x = []
        all_y = []
        all_parents = []
        for cat, l in d.iteritems():
            parents = l.keys()
            all_parents += parents

            expr_vals = l.values()
            x = [k[0] for k in expr_vals]
            y = [k[1] for k in expr_vals]
            all_x += x
            all_y += y

            ax_dict[rel_dir].scatter(
                x, y,
                s=10, color=color_dict[cat], lw=0, alpha=0.5,
                label=label_dict[cat]
            )
            ax_dict[rel_dir].set_xlim((-0.5,5))
            ax_dict[rel_dir].set_ylim(-0.5,5)
            ax_dict[rel_dir].set_aspect('equal')

        if report_outliers:
            upper_c, lower_c, outliers = find_outliers(all_x, all_y, all_parents)
            for ol in outliers:
                print(ol)
            ax_dict[rel_dir].plot([0,5], [upper_c, 5+upper_c], 'k--')
            ax_dict[rel_dir].plot([0,5], [lower_c, 5+lower_c], 'k--')

    if log_tpm:
        axarr[-1].set_xlabel("B cell parent log10(TPM + 1)", size='small')
        for ax in axarr:
            ax.set_ylabel("T cell parent log10(TPM + 1)", size='small')
    else:
        axarr[-1].set_xlabel("B cell parent TPM", size='small')
        for ax in axarr:
            ax.set_ylabel("T cell parent TPM", size='small')

    axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize='small')


    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("B_spec_mbed")
    parser.add_argument("T_spec_mbed")
    parser.add_argument("parent_expr_db")
    parser.add_argument("strand_combo_json")
    parser.add_argument("intrachrom_only", type=int)
    parser.add_argument("log_tpm", type=int)
    parser.add_argument("report_outliers", type=int)
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    with open(args.strand_combo_json, 'r') as f:
        strand_combo_dict = json.load(f)

    plot_data = get_plot_data(
        args.B_spec_mbed,
        args.T_spec_mbed,
        args.parent_expr_db,
        strand_combo_dict,
        args.intrachrom_only,
        args.log_tpm
    )

    draw_scatter(plot_data, args.log_tpm, args.report_outliers, args.plot_title, args.img_fname)
