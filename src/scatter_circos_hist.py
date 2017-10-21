'''
Take two circos hist data files with identical regions and do a scatter plot of
the height values, one point per region.
'''

import subprocess
import matplotlib.pyplot as plt

def check_regions_identical(fnames):
    region_sets = []
    cmd_str = "awk '{print $1,$2,$3}' %s | sort | uniq"
    for fn in fnames:
        rs = tuple(set(subprocess.check_output(cmd_str % fn, shell=True).strip().split("\n")))
        region_sets.append(rs)
    region_sets = tuple(set(region_sets))
    assert len(region_sets)==1, "Input files do not have identical regions!"

def get_plot_data(input_fnames, input_labels):
    input_fnames = input_fnames.split(",")
    input_labels = input_labels.split(",")

    check_regions_identical(input_fnames)

    plot_data = {l:{} for l in input_labels}
    for fn, label in zip(input_fnames, input_labels):
        with open(fn, 'r') as f:
            for line in f:
                line_list = line.strip().split()
                region = ":".join(line_list[0:3])
                data = float(line_list[3])
                plot_data[label][region] = data

    return plot_data

def draw_scatter(plot_data, input_labels, plot_title, img_fname):
    fig = plt.figure()
    plt.suptitle(plot_title)
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    input_labels = input_labels.split(",")

    regions = sorted(plot_data[input_labels[0]].keys())
    x = [plot_data[input_labels[0]][r] for r in regions]
    y = [plot_data[input_labels[1]][r] for r in regions]

    ax.scatter(x, y, c='red', s=4, lw=0)
    ax.set_xlabel(input_labels[0])
    ax.set_ylabel(input_labels[1])

    plt.savefig(img_fname)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fnames")
    parser.add_argument("input_labels")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(args.input_fnames, args.input_labels)

    draw_scatter(plot_data, args.input_labels, args.plot_title, args.img_fname)
