'''
Draw histograms of methylation level, mean, and variance for each set of regions
in the given data.
'''

import sys
import json
from itertools import product
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.stats import gaussian_kde

def get_plot_data(meth_data_json):
    with open(meth_data_json, 'r') as f:
        meth_data = json.load(f)

    data_types = ["level", "mean", "variance"]
    region_sets = meth_data.keys()
    meth_sets = list(
        set(
            reduce(lambda a,x: a+x, [v.keys() for v in meth_data.values()], [])
        )
    )

    plot_keys = ["/".join(t) for t in product(region_sets, meth_sets)]

    plot_data = {dt:{pk:[] for pk in plot_keys} for dt in data_types}

    for rs, d in meth_data.iteritems():
        for ms, e in d.iteritems():
            plot_key = rs + "/" + ms
            for region, f in e.iteritems():
                if f is None:
                    continue
                for dt, val in f.iteritems():
                    plot_data[dt][plot_key].append(val)

    return plot_data

def draw_histograms(plot_data, plot_title, img_fname):

    data_types = ["level", "mean", "variance"]

    fig, axarr = plt.subplots(len(data_types), 1, sharey=True)
    fig.set_size_inches(4, 15)

    ax_dict = {dt:ax for dt,ax in zip(data_types, axarr)}

    for dt, d in plot_data.iteritems():
        xmax = max(reduce(lambda a,x: a+x, d.values(), []))
        xvals = linspace(0, xmax, 1000)
        for pk, l in d.iteritems():
            density = gaussian_kde(l)
            yvals = density(xvals)
            ax_dict[dt].plot(xvals, yvals, label=pk)

    axarr[0].legend()

    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("meth_data_json")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(args.meth_data_json)

    draw_histograms(plot_data, args.plot_title, args.img_fname)
