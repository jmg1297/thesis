#Script that takes an "mBED" file (i.e., the output of bedtools intersect or bedtools closest; essentially two BED files next to each other with an extra column) and produces a heatmap showing the
#amount of each kind of feature in the second BED file there is in each of the features in the first BED file

import sys
import re
import subprocess
import argparse
from math import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib import colors
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

def get_headers(mbed_fname):
	'''
	Get the different categories of features the transcripts were intersected
	with
	'''
	cmd_str = "awk '{print $10}' %s | cut -d: -f1,2 | sort | uniq | sed '/\?/d'" % mbed_fname
	headers = sorted([h.strip() for h in subprocess.check_output(cmd_str, shell = True).strip().split("\n")])
	return headers

def get_props(raw_dict, headers):
	out_dict = {h:float(raw_dict[h])/raw_dict["length"] for h in headers}
	out_dict["UNIQ"] = max(0, 1 - sum(out_dict.values()))
	out_dict["tpm"] = np.log10(raw_dict["tpm"])
	out_dict["length"] = np.log10(raw_dict["length"])
	out_dict["coords"] = raw_dict["coords"]
	return out_dict

def mbed_to_dict(mbed_fname, headers):
	'''
	Parse the mBED file and create a dictionary containing the data we need for
	the heatmap, with the following structure:
	KEY transcript_id:
		-> KEY total_transcript_length : VALUE
		-> KEY rep_class : VALUE num_bases
		-> KEY rep class : VALUE num_bases
		...
	'''
	print("Parsing intersection file")
	sys.stdout.flush()

	num_lines = int(
		subprocess.check_output(
			"wc -l %s | awk '{print $1}'" % mbed_fname,
			shell = True
		)
	)
	percent = 0.0
	print('['+' '*50+']'),
	print('\r'),
	sys.stdout.flush()

	with open(mbed_fname, 'r') as mbed:
		raw_dict = {}
		for i,line in enumerate(mbed):

			percent = 100.0*(float(i)/num_lines)
			num_bars = int(floor(percent/2))
			print('\r'),
			print('['+'|'*num_bars+' '*(50-num_bars)+']'),
			sys.stdout.flush()

			line_list = re.split("\s+", line.strip())

			transcript_id = line_list[3]
			feature = re.sub("\?", "", ":".join(line_list[9].split(":")[0:2]))
			overlap = int(line_list[-1])

			try:
				raw_dict[transcript_id][feature] += overlap
			except KeyError:
				tpm = float(line_list[4])
				chrom = line_list[0]
				start = int(line_list[1])
				end = int(line_list[2])
				length = end - start

				raw_dict[transcript_id] = {h:0 for h in headers}
				raw_dict[transcript_id].update({
					"tpm":tpm,
					"length":length,
					"coords":[chrom, start, end]
				})
				try:
					raw_dict[transcript_id][feature] += overlap
				except KeyError:
					continue

	print('\r'),
	print('[' + '|'*50 + ']')
	print("Done")
	sys.stdout.flush()

	print("Creating dict with proportions")
	sys.stdout.flush()

	prop_dict = {k:get_props(raw_dict[k], headers) for k in raw_dict.keys()}

	print("Done")
	sys.stdout.flush()
	return prop_dict

def get_data_arrays(data_dict, headers):
	'''
	Convert the dict of data indexed by transcript into a set of arrays
	containing all of the data
	'''
	labels = headers + ["UNIQ"]
	transcripts = sorted(data_dict.keys())

	first = transcripts[0]
	d = data_dict[first]
	tpm_list = [np.log10(d["tpm"])]
	length_list = [np.log10(d["length"])]
	coords_list = [d["coords"]]
	feat_array = np.array([d[label] for label in labels])

	for k in transcripts[1:]:
		d = data_dict[k]
		tpm_list.append(d["tpm"])
		length_list.append(d["length"])
		coords_list.append(d["coords"])
		feat_array = np.vstack((feat_array, np.array([d[label] for label in labels])))

	tpm_array = np.array(tpm_list)
	length_array = np.array(length_list)

	return tpm_array, length_array, coords_list, feat_array

def make_color_list(headers):
	headers_dict = {}
	for h in headers:
		kv_pair = h.split(":")
		try:
			headers_dict[kv_pair[0]].append(kv_pair[1])
		except KeyError:
			headers_dict[kv_pair[0]] = [kv_pair[1]]
	families = sorted(headers_dict.keys())
	# TODO: Generalise away from retrotransposons
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

def cmp_chrs(x,y):
	#Custom comparison function for sorting chromosomes. Puts them in order 1,2,...,19, X, Y
	try:
		x = int(x)
	except ValueError:
		x = ord(x)
	try:
		y = int(y)
	except ValueError:
		y = ord(y)
	if x < y:
		return -1
	elif x == y:
		return 0
	elif x > y:
		return 1

def cmp_coords(x,y):
	# x and y are 3-element lists [chromosome, start, end]
	# Chromosomes have the order 1 - 19, X, Y
	# Coords are in numerical order
	# Uses the cmp_chrs function to order by chromosome first. In the case of a
	# tie (i.e., same chromosome), starting coordinate is used to break the tie
	# (unless they start at the same place).
	chr_cmp = cmp_chrs(x[0], y[0])

	if chr_cmp < 0:
		return -1
	elif chr_cmp == 0:
		if x[1] < y[1]:
			return -1
		elif x[1] == y[1]:
			return 0
		elif x[1] > y[1]:
			return 1
	elif chr_cmp > 0:
		return 1

def get_coord_order(coords_list):
	# Returns the indices of the coordinates ion coordList in sorted order.
	# Sort uses a custom function for comparing coordinates
	return [i for i,c in sorted(enumerate(coords_list), key = lambda x: x[1], cmp = cmp_coords)]

def make_heatmap(feat_array, tpm_array, length_array, coords_list, headers, sort_by, bar_right, fig_fname):
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

	num_transcripts = len(coords_list)

	print("Sorting data by {}".format(sort_by))
	sys.stdout.flush()

	if sort_by == "cluster":
		# Do the clustering
		Z = linkage(feat_array, 'ward')

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

		# Sort the data according to the leaves list
		feat_array = feat_array[leaves_list(Z)]
		tpm_array = tpm_array[leaves_list(Z)]
		length_array = length_array[leaves_list(Z)]

	elif sort_by == "chromosome":
		# Order the data according to chromosome and starting coordinate
		coord_order = get_coord_order(coords_list)

		feat_array = feat_array[coord_order]
		tpm_array = tpm_array[coord_order]
		length_array = length_array[coord_order]

		# Put chromosome numbers to the left of the heatmap
		chrom_axis = fig.add_axes([
			0.08/tot_width,
			0.08/tot_height,
			1.8/tot_width,
			8/tot_height
		])

		all_chroms = [c[0] for c in coords_list]
		chroms = sorted(list(set(all_chroms)), cmp = cmp_chrs)
		chrom_axis.set_ylim([0,1])
		chrom_axis.set_xlim([0,1])
		chrom_props = {chrom:float(all_chroms.count(chrom))/len(all_chroms) for chrom in chroms}
		y = 0
		for chrom in chroms:
			cmap = matplotlib.cm.get_cmap("Purples")
			col = cmap(float(chrom_props[chrom])/max(chrom_props.values()))
			chrom_axis.add_patch(patches.Rectangle(
				(0,y),
				1,
				chrom_props[chrom],
				facecolor = col
			))
			plt.text(
				x = 0.5,
				y = y + float(chrom_props[chrom])/2,
				s = str(chrom),
				verticalalignment = "center",
				horizontalalignment = "center"
			)
			y += chrom_props[chrom]

	elif sort_by == "tpm":

		#Order the data according to the TPM
		order = np.argsort(tpm_array)
		feat_array = feat_array[list(order)]
		length_array = length_array[order]
		tpm_array = tpm_array[order]

		# Put a horizontal barchart of TPMs to the left of the heatmap
		tpm_axis = fig.add_axes([
			0.08/tot_width,
			0.08/tot_height,
			1.8/tot_width,
			8/tot_height
		])
		tpm_axis.barh(
			range(1, num_transcripts+1),
			tpm_array,
			height = 1,
			align = "edge",
			lw = 0
		)
		tpm_axis.set_ylim([1, num_transcripts+1])
		tpm_axis.set_yticklabels([])
		tpm_axis.tick_params(axis = "x", labelsize = 10)
		tpm_axis.set_xlabel("log10(TPM)")

	elif sort_by == "length":
		#Order the data according to their lengths
		order = np.argsort(length_array)
		feat_array = feat_array[order]
		tpm_array = tpm_array[order]
		length_array = np.sort(length_array)

		#Put a horizontal barchart of lengths to the left of the heatmap
		length_axis = fig.add_axes([
			0.08/tot_width,
			0.08/tot_height,
			1.8/tot_width,
			8/tot_height
		])
		length_axis.barh(
			range(1, num_transcripts+1),
			length_array,
			height = 1,
			align = "edge",
			lw = 0
		)
		length_axis.set_ylim([1, num_transcripts+1])
		length_axis.set_yticklabels([])
		length_axis.tick_params(axis = "x", labelsize = 10)
		length_axis.set_xlabel("log10(Length)")

	print("Finished sorting")
	# Convert the raw data into a matrix that can be used by imshow.
	# Each transcript is divided into 100 bins, and each class of repeat is
	# assigned a number of bins corresponding to the proportion of the sequence
	# it contributes

	print("Creating heatmap matrix")
	sys.stdout.flush()

	print('['+' '*50+']'),
	print('\r'),
	percent = 0.0

	num_classes = len(headers) + 1
	is_first_row = True
	for i,row in enumerate(feat_array):

		percent = 100.0*(float(i)/num_transcripts)
		num_bars = int(floor(percent/2))
		print('\r'),
		print('['+'|'*num_bars+' '*(50-num_bars)+']'),
		sys.stdout.flush()


		if -1 in np.sign(row):
			print("removing neg row {}".format(i))
			print(row)
			continue
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
		1.96/tot_width,
		0.08/tot_height,
		8/tot_width,
		8/tot_height
	])
	hmap_axis.set_xlabel("% Transcript Length")

	# Use imshow to create the heatmap from the plot_data array
	img = hmap_axis.imshow(
		plot_data,
		aspect = "auto",
		interpolation = "nearest",
		origin = "lower",
		cmap = cmap,
		norm = norm
	)

	# Remove the y axis labels
	img.axes.get_yaxis().set_visible(False)
	labels = headers + ["UNIQ"]

	# Add the color bar
	cbar_axis = fig.add_axes([
		1.96/tot_width,
		8.31/tot_height,
		8/tot_width,
		0.25/tot_height
	])
	cbar = plt.colorbar(
		img,
		cax = cbar_axis,
		cmap = cmap,
		norm = norm,
		boundaries = bounds,
		ticks = range(len(labels)),
		orientation = "horizontal"
	)
	cbar_labels = [l.split(":")[1] for l in labels[0:-1]] + ["UNIQ"]
	cbar.ax.set_xticklabels(cbar_labels, fontsize = 5)

	# Add the horizontal barplot of TPMs or lengths
	if bar_right != "none":
		bar_axis = fig.add_axes([
			10.04/tot_width,
			0.08/tot_height,
			2/tot_width,
			8/tot_height
		])

		if bar_right == 'tpm':
			bar_axis.barh(
				range(1, num_transcripts+1),
				tpm_array,
				height = 1,
				align = "edge",
				lw = 0
			)
			bar_axis.set_xlabel("log10(TPM)")

		elif bar_right == 'length':
			bar_axis.barh(
				range(1, num_transcripts+1),
				length_array,
				height = 1,
				align = "edge",
				lw = 0
			)
			bar_axis.set_xlabel("log10(Length)")

		bar_axis.set_ylim([1, num_transcripts+1])
		bar_axis.set_yticklabels([])
		bar_axis.tick_params(axis = "x", labelsize = 10)

	# Save the figure to the file name given
	print("Done, saving figure to {}".format(fig_fname))
	sys.stdout.flush()
	plt.savefig(fig_fname, bbox_inches = "tight")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description = "Create a heatmap of transcript sequence content from \
						the output of bedtools intersect"
	)
	parser.add_argument(
		"--intersectFilename",
		help = "File of intersections between transcripts and repeats (or other features)",
		dest = "intersect_fname",
		type = str,
		required = True
	)
	parser.add_argument(
		"--outputImg", help = "Name of the file to save the heatmap to",
		dest = "fig_fname",
		type  = str,
		required = True
	)
	parser.add_argument(
		"--sortBy",
		help = """
			Determines the order in which transcripts are arranged in the heatmap.
			Options are \"cluster\" (default), \"tpm\", \"length\", \"chromosome\".
			Selecting cluster will display a dendrogram on the left of the figure.
			Selecting chromosome will display chromosome numbers to the left of
			the heatmap.
			Selecting tpm or length will display a horizontal barchart of the
			respective variables to the left of the heatmap. A separate option
			controls the barchart to the right of the heatmap.
		""",
		dest = "sort_by",
		default = "cluster",
		type = str,
		choices = ["cluster", "tpm", "length", "chromosome"]
	)
	parser.add_argument(
		"--barRight",#
		help = """
			Controls which variable (if any) is plotted to the right of the
			heatmap.
			Options are \"none\" (default), \"tpm\", \"length\". Note the latter
			two options will be displayed after being transformed to log10 values
		""",
		dest = "bar_right",
		default = "none",
		type = str,
		choices = ["none", "length", "tpm"]
	)

	args = parser.parse_args()

	headers = get_headers(args.intersect_fname)

	data_dict = mbed_to_dict(args.intersect_fname, headers)

	tpm_array, length_array, coords_list, feat_array = get_data_arrays(data_dict, headers)

	make_heatmap(
		feat_array,
		tpm_array,
		length_array,
		coords_list,
		headers,
		args.sort_by,
		args.bar_right,
		args.fig_fname
	)
