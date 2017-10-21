'''
Script to produce a stacked horizontal barchart showing contributions of
different RT families to sets of regions. In particular, used to show which RTs
have gaps, expressed gaps, etc.

Assume that all BED file names describe a rmsk RT element in the colon-separated
format used in this project, e.g., LINE:L1:L1_Mur3:10. Extra fields after this
are allowed.
'''

import sys
import re
import subprocess
import logging

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

logging.basicConfig(level=logging.INFO)

def get_proportions(bed_fname):
	logging.info("Getting proportions from {} ... ".format(bed_fname)),

	raw_dict = {"LINE":{}, "SINE":{}, "LTR":{}}
	prop_dict = {}

	total = int(
		subprocess.check_output(
			"wc -l %s | awk '{print $1}'" % bed_fname,
			shell=True
		)
	)

	with open(bed_fname, 'r') as f:
		for line in f:
			line_list = line.strip().split()
			name_list = line_list[3].split(":")
			rt_type = name_list[0]
			family = name_list[1]
			try:
				raw_dict[rt_type][family] += 1
			except KeyError:
				raw_dict[rt_type][family] = 1

	for rt_type in raw_dict.keys():
		prop_dict[rt_type] = {
			family:100*float(raw_dict[rt_type][family])/total
				for family in raw_dict[rt_type].keys()
		}

	logging.info("Done")

	return prop_dict

def get_prop_dicts(bed_fnames, labels):
	out_dict = {}
	for fn,l in zip(bed_fnames, labels):
		out_dict[l] = get_proportions(fn)
	return out_dict

def get_color_dict(rt_types_dict):
	'''
	Assume for simplicity we are dealing only with families of SINEs, LINEs,
	and LTRs
	'''

	cmap_dict = {
		"LINE":matplotlib.cm.get_cmap("Greens"),
		"SINE":matplotlib.cm.get_cmap("Blues"),
		"LTR":matplotlib.cm.get_cmap("Reds")
	}

	color_dict = {"LINE":{}, "SINE":{}, "LTR":{}}

	for rt_type in rt_types_dict.keys():
		cmap = cmap_dict[rt_type]
		num_fams = len(rt_types_dict[rt_type])
		counter = 1.0
		for family in sorted(rt_types_dict[rt_type]):
			color_dict[rt_type][family] = cmap(counter/num_fams)
			counter += 1
	return color_dict

def draw_prop_barcharts(prop_dicts, labels, plot_title, img_fname):
	'''
	Input should be a dictionary where keys are file labels and values are
	dictionaries of proportions
	'''

	num_files = len(labels)

	fig = plt.figure(figsize=[12, 3*num_files])
	plt.suptitle(plot_title)

	axis = fig.add_axes([0.05,0.05,0.9,0.9])

	axis.set_ylabel("Region Set")
	axis.set_xlabel("Percentage")

	height = 0.8

	axis.set_xlim([-5,105])
	axis.set_ylim([0.5, num_files+1])
	axis.set_yticks([i + height/2 for i in range(1, num_files+1)])
	axis.set_yticklabels(labels)

	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)

	rt_types_dict = {"LINE":[], "SINE":[], "LTR":[]}
	for l in labels:
		for rt_type in prop_dicts[l].keys():
			for family in prop_dicts[l][rt_type].keys():
				rt_types_dict[rt_type].append(family)

	rt_types_dict = {
		k:sorted(list(set(rt_types_dict[k]))) for k in rt_types_dict.keys()
	}

	color_dict = get_color_dict(rt_types_dict)

	for i, label in enumerate(labels):
		prop_dict = prop_dicts[label]
		sum_width = 0
		for rt_type in sorted(prop_dict.keys()):
			for family in prop_dict[rt_type].keys():
				percentage = prop_dict[rt_type][family]
				axis.barh(
					i+1,
					percentage,
					height=height,
					left=sum_width,
					color=color_dict[rt_type][family],
					lw=0
				)

				if percentage > 5:
					axis.text(
						sum_width + percentage/2,
						i + 1 + height/2,
						"%s:%s\n%.2f" % (rt_type, family, percentage),
						horizontalalignment='center',
						verticalalignment='center',
						fontsize=10
					)
				sum_width += percentage


	legend_patches = []
	for rt_type in sorted(color_dict.keys()):
		legend_patches \
			+= [
				mpatches.Patch(
					color=color_dict[rt_type][family],
					label=rt_type+":"+family
				) for family in sorted(color_dict[rt_type].keys()
			)]

	plt.legend(
		handles=legend_patches,
		loc=2,
		bbox_to_anchor=(1.05, 1.05),
		prop={'size':10}
	)

	plt.savefig(img_fname, bbox_inches='tight')

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--bedFname", action="append")
	parser.add_argument("-l", "--label", action="append")
	parser.add_argument("plot_title")
	parser.add_argument("img_fname")
	args = parser.parse_args()

prop_dicts = get_prop_dicts(args.bedFname, args.label)

draw_prop_barcharts(prop_dicts, args.label, args.plot_title, args.img_fname)
