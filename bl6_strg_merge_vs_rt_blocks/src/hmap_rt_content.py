'''
Given a BED file of transcripts with summarised retrotransposon overlaps,
draw a heatmap showing the content of each transcript, clustered based on
content.
'''

import subprocess
import sys
import logging
import re
from math import *

import numpy as np
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib import colors
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

logging.basicConfig(level=logging.INFO)

def get_headers(summary_bed_fname):
    '''
    Get the retrotransposon types present in the intersection
    '''
    cmd_str = "cat %s | awk '{print $7}' | tr ';' '\n' | cut -d, -f1 | sed -r 's/\?//g' | sort | uniq" % summary_bed_fname
    headers = [
        h.strip() for h in subprocess.check_output(
            cmd_str, shell=True).strip().split("\n")]
    return headers

def get_exon_length_dict(exon_bed_fname):
    exon_length_dict = {}
    with open(exon_bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = ".".join(line_list[3].split(".")[0:3])
            length = int(line_list[2]) - int(line_list[1])
            try:
                exon_length_dict[transcript] += length
            except KeyError:
                exon_length_dict[transcript] = length
    return exon_length_dict

def get_intersect_dict(summary_bed_fname, exon_bed_fname, headers):
    count_dict = {}
    length_dict = get_exon_length_dict(exon_bed_fname)
    labels = headers + ["UNIQ"]

    with open(summary_bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = line_list[3]
            count_dict[transcript] = {l:0 for l in labels}
            overlaps = line_list[6].split(";")
            for overlap in overlaps:
                label, bp = overlap.split(",")
                label = re.sub("\?", "", label)
                count_dict[transcript][label] += int(bp)
            count_dict[transcript]["UNIQ"] = max(0, length_dict[transcript] - sum(count_dict[transcript].values()))

    prop_dict = {
        transcript:{
            l:float(data[l])/length_dict[transcript]
            for l in labels}
        for transcript, data in count_dict.iteritems()}

    return prop_dict

def get_data_array(summary_bed_fname, exon_bed_fname, headers):

    intersect_dict = get_intersect_dict(summary_bed_fname, exon_bed_fname, headers)

    labels = headers + ["UNIQ"]
    first = True

    for transcript, data in intersect_dict.iteritems():
        if first:
            feat_array = np.array([data[l] for l in labels])
            first = False
        else:
            feat_array = np.vstack((feat_array, np.array([data[l] for l in labels])))

    num_transcripts = len(intersect_dict)

    return feat_array, num_transcripts

def make_color_list(headers):
	headers_dict = {}
	for h in headers:
		kv_pair = h.split(":")
		try:
			headers_dict[kv_pair[0]].append(kv_pair[1])
		except KeyError:
			headers_dict[kv_pair[0]] = [kv_pair[1]]
	families = sorted(headers_dict.keys())
	color_list = []

	for i,f in enumerate(families):
		sub_fams = sorted(headers_dict[f])
		num_colors = len(sub_fams)
		if f == "LINE":
			cmap = matplotlib.cm.get_cmap("Greens")
		elif f == "SINE":
			cmap = matplotlib.cm.get_cmap("Blues")
		elif f == "LTR":
			cmap = matplotlib.cm.get_cmap("Reds")
		for j in range(num_colors):
			color_list.append(cmap(float(j+1)/num_colors))
	color_list.append("white")
	return color_list

def make_heatmap(feat_array, num_transcripts, headers, img_fname):
    '''
    Draw a heatmap showing the relative amount of each type of feature in each
    transcript.

    Figure layout:
    0.08 inch borders around edges
    0.08 inch gaps between different plots, except between top of heatmap and the colorbar
    [left edge, bottom edge, width, height]
    Dendrogram: 0.08, 0.08, 1.8, 8
    Heatmap: 1.96, 0.08, 8, 8
    Barchart: 10.04, 0.08, 2, 8
    Cbar: 1.96, 8.31, 8, 0.25

    Total Width = 12.12
    Total Height = 8.64
    '''
    # TODO: automate these calculations
    tot_width = 12.12
    tot_height = 8.64

    # Create the figure object that will be populated and saved to a file
    fig = plt.figure(figsize = (tot_width, tot_height), dpi = 100)

    logging.info("Clustering data")

    # Do the clustering
    Z = linkage(feat_array, 'ward')

    logging.info("Finished clustering")

    # Sort the data according to the leaves list
    #feat_array = feat_array[leaves_list(Z)]

    #logging.info("Finished sorting")

    logging.info("Plotting dendogram")

    # Create the dendrogram
    # Add an axis for the dendrogram
    dendro_axis = fig.add_axes([
        0.08/tot_width,
        0.08/tot_height,
        1.8/tot_width,
        8/tot_height
    ])

    # Plot the dendrogram
    dendrogram(Z, orientation = "right")

    # Remove x and y tick labels and add a y axis label
    dendro_axis.set_xticks([])
    dendro_axis.set_yticks([])
    dendro_axis.set_ylabel("Transcripts (%d)" % num_transcripts)

    # Convert the raw data into a matrix that can be used by imshow.
    # Each transcript is divided into 100 bins, and each class of repeat is
    # assigned a number of bins corresponding to the proportion of the sequence
    # it contributes

    logging.info("Creating heatmap matrix")

    print('['+' '*50+']'),
    print('\r'),
    percent = 0.0

    num_classes = len(headers) + 1
    is_first_row = True
    #for i,row in enumerate(feat_array):
    for c, idx in enumerate(leaves_list(Z)):

        row = feat_array[idx]

        percent = 100.0*(float(c)/num_transcripts)
        num_bars = int(floor(percent/2))
        print('\r'),
        print('['+'|'*num_bars+' '*(50-num_bars)+']'),
        sys.stdout.flush()

        if -1 in np.sign(row):
            print("Neg row {}".format(idx))
            print(row)
            sys.exit()
        # Create the new row by iteratively extending a list with the
        # appropriate classification repeated however many times
        new_row = []
        for j in range(num_classes):
            new_row += [j for i in range(int(round(100*row[j])))]

        # If rounding errors have created too many or too few columns,
        # remove or add as appropriate
        if len(new_row) > 100:
            new_row = new_row[0:100]
        elif len(new_row) < 100:
            if len(new_row) >= 95:
                while len(new_row) != 100:
                    new_row = np.append(new_row, new_row[-1])
            else:
                print("Insufficient columns")
                print(row)
                print(new_row)
                sys.exit()
                continue

        # If this is the first row, create the new numpy array; otherwise
        # concatenate the new row to the existing array
        if is_first_row:
            plot_data = np.array(new_row)
            is_first_row = False
        else:
            plot_data = np.vstack((plot_data, new_row))

    print('\r'),
    print('['+'|'*50+']')
    print("Done")
    print("Drawing heatmap")
    sys.stdout.flush()

    # Create the colormap
    cmap = colors.ListedColormap(make_color_list(headers))

    bounds = list(np.arange(-0.5, num_classes+0.5, 1))
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Add an axis to the figure for the heatmap
    hmap_axis = fig.add_axes([
        1.96/tot_width,    0.08/tot_height, 8/tot_width, 8/tot_height])

    hmap_axis.set_xlabel("% Transcript Length")

    # Use imshow to create the heatmap from the plot_data array
    img = hmap_axis.imshow(
        plot_data,
        aspect="auto", interpolation="nearest", origin="lower",
        cmap=cmap, norm=norm)

    # Remove the y axis labels
    img.axes.get_yaxis().set_visible(False)
    labels = headers + ["UNIQ"]

    # Add the color bar
    cbar_axis = fig.add_axes([
        1.96/tot_width, 8.31/tot_height, 8/tot_width, 0.25/tot_height])

    cbar = plt.colorbar(
        img,
        cax=cbar_axis, cmap=cmap, norm=norm,
        boundaries=bounds, ticks=range(len(labels)),
        orientation="horizontal")

    cbar_labels = [l.split(":")[1] for l in labels[0:-1]] + ["UNIQ"]
    cbar.ax.set_xticklabels(cbar_labels, fontsize=12, rotation=60, ha='left')
    cbar_axis.xaxis.tick_top()
    cbar_axis.xaxis.set_label_position('top')

    # Save the figure to the file name given
    logging.info("Finished heatmap, saving figure to {}".format(img_fname))
    plt.savefig(img_fname, bbox_inches="tight")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("summary_bed_fname")
    parser.add_argument("exon_bed_fname")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    headers = get_headers(args.summary_bed_fname)

    feat_array, num_transcripts = get_data_array(args.summary_bed_fname, args.exon_bed_fname, headers)

    make_heatmap(feat_array, num_transcripts, headers, args.img_fname)
