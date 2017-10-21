'''
Plot barcharts showing the relative proportions of each chromosome for expressed
vs all retrocopies, and same for retrocopy parents
'''

import subprocess
import tempfile as tf
import matplotlib.pyplot as plt
import numpy as np
import logging

logging.basicConfig(level=logging.DEBUG)

def get_chrom_props(bed_fname):
    cmd_str = "cat %s | awk '{print $1}' | sort | uniq -c | awk '{print $2,$1}'" % bed_fname
    chrom_counts = subprocess.check_output(cmd_str, shell=True).strip().split("\n")
    total = sum([int(l.split()[1]) for l in chrom_counts])
    return {chrom:float(count)/total for chrom,count in map(lambda x: x.split(), chrom_counts)}

def get_stats(data):
    return {"mean":np.mean(data), "std":np.std(data)}

def get_random_data(bed_fname, sample_size, num_samples):
    chroms = subprocess.check_output(
        "cat %s | awk '{print $1}' | sort | uniq" % bed_fname,
        shell=True
    ).strip().split("\n")

    data = {c:[] for c in chroms}

    for i in range(num_samples):
        if i%100 == 0:
            logging.info("Sample {}/{}".format(i, num_samples))

        with tf.NamedTemporaryFile() as tmp:
            subprocess.call("shuf -n %d %s > %s" % (sample_size, bed_fname, tmp.name), shell=True)
            cp = get_chrom_props(tmp.name)
            for c,p in cp.iteritems():
                data[c].append(p)

    logging.info("Got all random samples")

    random_data = {c:get_stats(data[c]) for c in chroms}

    return random_data

def get_plot_data(expr_retrocopies_bed, all_retrocopies_bed, expr_rc_parents_bed, all_rc_parents_bed, sample_size, num_samples):

    plot_data = {
        "retrocopies":{
            "expressed":get_chrom_props(expr_retrocopies_bed),
            "all":get_random_data(all_retrocopies_bed, sample_size, num_samples)
        },
        "parents":{
            "expressed":get_chrom_props(expr_rc_parents_bed),
            "all":get_random_data(all_rc_parents_bed, sample_size, num_samples)
        }
    }

    return plot_data

def get_chroms(*args):
    return subprocess.check_output(
        "cat %s | awk '{print $1}' | sort | uniq" % " ".join(args),
        shell=True).strip().split("\n")

def cmp_chrs(x, y):
    if "chr" in x:
        x = x[3:]

    if "chr" in y:
        y = y[3:]

    try:
        x = int(x)
    except ValueError:
        x = ord(x)

    try:
        y = int(y)
    except ValueError:
        y = ord(y)

    return x-y

def draw_barcharts(plot_data, chroms, plot_title, img_fname):
    width = 0.4
    color_dict = {"expressed":"red", "all":"gray"}
    offset_dict = {"expressed":0, "all":width}

    chroms = sorted(chroms, cmp=cmp_chrs)

    fig, axarr = plt.subplots(2, 1, sharey = True)
    fig.set_size_inches(10, 10)

    ax_dict = {"retrocopies":axarr[0], "parents":axarr[1]}

    ind = np.arange(len(chroms))

    for plot, pdata in plot_data.iteritems():
        for status, chrom_data in pdata.iteritems():
            heights = []
            errors = []
            if status == "expressed":
                for c in chroms:
                    try:
                        heights.append(chrom_data[c])
                    except KeyError:
                        heights.append(0.0)
                ax_dict[plot].bar(ind + offset_dict[status], heights, width, color=color_dict[status], alpha=0.7, lw=0)
            else:
                for c in chroms:
                    try:
                        heights.append(chrom_data[c]["mean"])
                        errors.append(chrom_data[c]["std"])
                    except KeyError:
                        heights.append(0.0)
                        errors.append(0.0)
                ax_dict[plot].bar(ind + offset_dict[status], heights, width, yerr=errors, color=color_dict[status], alpha=0.7, lw=0)
    '''

    heights1 = []
    errors1 = []
    for c,d in plot_data["retrocopies"]["all"].iteritems():
        heights1.append(d["mean"])
        errors1.append(d["std"])

    axarr[0].bar(ind, heights1, width, yerr=errors1, color="grey")

    heights2 = []
    errors2 = []
    for c,d in plot_data["parents"]["all"].iteritems():
        heights2.append(d["mean"])
        errors2.append(d["std"])

    axarr[0].bar(ind + width, heights2, width, yerr=errors2, color="red")

    '''
    for ax in axarr:
        ax.set_ylim((0,0.15))
        ax.set_xticks(ind + width)
        ax.set_xticklabels(chroms, rotation=90)
        ax.set_xlim([-1, len(chroms)+1])
        ax.set_ylabel("Proportion")

    axarr[-1].set_xlabel("Chromosome")

    axarr[0].text(0, 0.9*max(axarr[0].get_ylim()), "Retrocopies", fontsize=16)
    axarr[1].text(0, 0.9*max(axarr[1].get_ylim()), "Parents", fontsize=16)

    plt.suptitle(plot_title)

    plt.savefig(img_fname)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("expr_retrocopies_bed")
    parser.add_argument("all_retrocopies_bed")
    parser.add_argument("expr_rc_parents_bed")
    parser.add_argument("all_rc_parents_bed")
    parser.add_argument("sample_size", type=int)
    parser.add_argument("num_samples", type=int)
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(
        args.expr_retrocopies_bed,
        args.all_retrocopies_bed,
        args.expr_rc_parents_bed,
        args.all_rc_parents_bed,
        args.sample_size,
        args.num_samples
    )

    chroms = get_chroms(args.all_retrocopies_bed, args.all_rc_parents_bed)

    draw_barcharts(plot_data, chroms, args.plot_title, args.img_fname)
