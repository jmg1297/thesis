'''
Given two TSV files of protein data with abundance values before and after
normalisation, draw density plots of the abundance values for each sample.
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

def get_protein_data(tsv_fname, abundance_fields, log=False):
    data = {f:[] for f in abundance_fields}
    with open(tsv_fname, 'r') as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            line_dict = {h:v for h,v in zip(headers, line.strip().split("\t"))}
            for f in abundance_fields:
                val = line_dict[f]
                try:
                    val = float(val)
                except ValueError:
                    continue

                if log:
                    val = np.log2(val)

                data[f].append(val)
    return data

def get_plot_data(raw_tsv_fname, normed_tsv_fname, abundance_fields):
    plot_data = {ds:{} for ds in ["raw", "normed"]}

    plot_data["raw"] = get_protein_data(
        raw_tsv_fname, abundance_fields, log=True)

    plot_data["normed"] = get_protein_data(normed_tsv_fname, abundance_fields)

    return plot_data

def draw_density_plots(plot_data, abundance_field, field_color, img_fname):
    color_dict = {f:c for f,c in zip(abundance_field, field_color)}

    fig, axarr = plt.subplots(1, 2, sharex=True, sharey=True)
    fig.set_size_inches(12,4)
    ax_dict = {"raw":axarr[0], "normed":axarr[1]}
    title_dict = {"raw":"Raw Values", "normed":"Normalised Values"}

    all_values = reduce(
        lambda a,x: a+x,
        [reduce(lambda a,x: a+x, d.values(), []) for d in plot_data.values()],
        [])

    xvals = np.linspace(min(all_values), max(all_values), 1000)

    for dataset, data in plot_data.iteritems():
        for field, values in data.iteritems():
            density = gaussian_kde(values)
            yvals = density(xvals)
            ax_dict[dataset].plot(
                xvals, yvals, color=color_dict[field], lw=1, label=field)
            ax_dict[dataset].set_title(title_dict[dataset])
            ax_dict[dataset].set_xlabel("$\log_2(abundance)$")
    axarr[0].set_ylabel("Density")

    plt.legend()

    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_tsv_fname")
    parser.add_argument("normed_tsv_fname")
    parser.add_argument("-f", "--abundance_field", action="append", type=str)
    parser.add_argument("-c", "--field_color", action="append", type=str)
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(
        args.raw_tsv_fname, args.normed_tsv_fname, args.abundance_field)

    draw_density_plots(
        plot_data, args.abundance_field, args.field_color, args.img_fname)
