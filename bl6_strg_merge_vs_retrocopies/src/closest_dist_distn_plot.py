'''
Plot the distribution of distance to the nearest upstream retrotransposon block
for each of several sets of regions. Also generate random regions using length
distribution from another file, find closest upstream RT blocks, and plot
results
'''

import sys
import re
import subprocess
import tempfile as tf
from numpy import linspace
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def generate_random_data(chrom_sizes_fname, length_distn_fname, feat_bed_fname):
    tmp_unsorted = tf.NamedTemporaryFile()
    tmp_bed = tf.NamedTemporaryFile()
    tmp_closest = tf.NamedTemporaryFile()

    cmd_str = "python $HOME/projects/src/get_random_regions.py {} {} 50 15000 {}"

    subprocess.call(
        cmd_str.format(chrom_sizes_fname, length_distn_fname, tmp_unsorted.name),
        shell=True
    )

    print("Sorting random regions")
    sys.stdout.flush()

    subprocess.call(
        "sort -k1,1 -k2,2n {} > {}".format(tmp_unsorted.name, tmp_bed.name),
        shell=True
    )

    print("Sorted")
    print("Getting closest RTs")
    sys.stdout.flush()

    cmd_str = "bedtools closest -d -a {} -b {} > {}"

    subprocess.call(
        cmd_str.format(tmp_bed.name, feat_bed_fname, tmp_closest.name),
        shell=True
    )

    print("Done")

    tmp_unsorted.close()
    tmp_bed.close()
    tmp_closest.seek(0)
    return tmp_closest

def get_dist(line):
    dist = int(line.strip().split()[-1])
    return abs(dist)

def get_plot_data(input_fnames, labels, rand_file):

    print("Getting plot data")

    assert "RANDOM" not in labels, "RANDOM is a reserved label"

    plot_data = {l:[] for l in labels + ["RANDOM"]}

    for fn, l in zip(input_fnames, labels):
        with open(fn, 'r') as f:
            for line in f:
                dist = get_dist(line)
                plot_data[l].append(dist)

    for line in rand_file:
        dist = get_dist(line)
        plot_data["RANDOM"].append(dist)

    rand_file.close()
    return plot_data

def draw_density_plot(plot_data, labels, colors, plot_title, img_fname):

    assert "black" not in colors, "Black is a reserved color"

    color_dict = {l:c for l,c in zip(labels, colors)}
    color_dict["RANDOM"] = "black"

    #xmax = max(reduce(lambda a,x: a+x, plot_data.values(), []))
    xmax = 200000 # hard-coded to improve plot - sorry future users
    xvals = linspace(0, xmax, 1000)

    fig, ax = plt.subplots(1,1)
    #fig.set_size_inches(5, 5)

    for label, data in plot_data.iteritems():
        density = gaussian_kde(data)
        density.covariance_factor = lambda : .1
        density._compute_covariance()
        yvals = density(xvals)
        ax.plot(xvals, yvals, color=color_dict[label], label=label)

    ax.set_ylabel("Density")
    ax.set_xlabel("Distance (bp)")

    ax.set_xlim([0,xmax])

    ax.legend()

    plt.suptitle(plot_title)

    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFname", action="append")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("-c", "--color", action="append")
    parser.add_argument("chrom_sizes_fname")
    parser.add_argument("length_distn_fname")
    parser.add_argument("feat_bed_fname")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    rand_file = generate_random_data(
        args.chrom_sizes_fname, args.length_distn_fname, args.feat_bed_fname
    )

    plot_data = get_plot_data(args.inputFname, args.label, rand_file)

    draw_density_plot(
        plot_data, args.label, args.color, args.plot_title, args.img_fname
    )
