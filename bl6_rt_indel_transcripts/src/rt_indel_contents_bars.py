# Script to visualise the number of RT gaps of each classification
# Show all RT gaps and just the expressed ones

import sys
import re
import os
import subprocess
from math import *
import json

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def get_proportions_dict(fnames, labels, classes):
	'''
	For each file of classified indels, get proportions of each classification
	'''

	counts_dict = {label:{c:0 for c in classes} for label in labels}
	proportions_dict = {label:{} for label in labels}
	totals = {label:0 for label in labels}

	for fn,label in zip(fnames, labels):

		print("Getting data from classification file {}".format(fn))
		print('[' + ' '*50 + ']'),
		print('\r'),
		num_lines = int(
			subprocess.check_output(
				"wc -l %s | awk '{print $1}'" % fn,
				shell=True
			)
		)
		percent = 0.0

		totals[label] = num_lines

		with open(fn, 'r') as f:
			for i,line in enumerate(f):
				line_list = re.split("\s+", line.strip())
				classification = line_list[6]
				try:
					counts_dict[label][classification] += 1
				except KeyError:
					if "RNA" in classification:
						counts_dict[label]["RNA"] += 1
					elif "?" in classification:
						classification = re.sub("\?", "", classification)
						counts_dict[label][classification] += 1
					else:
						counts_dict[label]["Other"] += 1

				percent = 100.0*(float(i)/num_lines)
				num_bars = int(floor(percent/2))
				print('\r'),
				print('[' + '|'*num_bars + ' '*(50-num_bars) + ']'),
				sys.stdout.flush()

			print('\r'),
			print('[' + '|'*50 + ']'),
			sys.stdout.flush()
			print("\nDone")

		proportions_dict[label] = {
			c:100*float(counts_dict[label][c])/num_lines for c in classes
		}

	return proportions_dict, totals

def draw_barchart(proportions_dict, totals, label_list, classes, plot_title, img_fname):
	'''
	Make a horizontal barchart where each bar represents one set of indels
	'''

	print("Drawing barchart ... "),
	sys.stdout.flush()

	color_dict = {
					"DNA":"brown",
					"LINE":"green",
					"Low_complexity":"cyan",
					"LTR":"red",
					"pseudogene":"purple",
					"RNA":"yellow",
					"Satellite":"magenta",
					"Simple_repeat":"orange",
					"SINE":"blue",
					"Other":"black",
					"Mixed":"grey",
					"NONE":"beige"
				}

	num_bars = len(proportions_dict)

	fig = plt.figure(figsize = [12,3*num_bars])
	plt.suptitle(plot_title)

	axis = fig.add_axes([0.05,0.05,0.9,0.9])

	axis.set_ylabel("Region set")
	axis.set_xlabel("Percentage of regions")

	height = 0.8

	axis.set_xlim([-5,105])
	axis.set_ylim([0.5, num_bars+1])
	axis.set_yticks([i + height/2 for i in range(1, num_bars+1)])
	ytick_labels = [
		"{} ({})".format(label, totals[label]) for label in label_list
	]
	axis.set_yticklabels(ytick_labels)

	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)

	for i,label in enumerate(label_list):
		data = proportions_dict[label]
		sum_width = 0
		for c in classes:
			percentage = data[c]
			axis.barh(
				i+1,
				percentage,
				height=height,
				left=sum_width,
				color=color_dict[c],
				lw=0,
				alpha=0.7
			)

			if percentage > 5:
				axis.text(
					sum_width + percentage/2,
					i + 1 + height/2,
					"%s\n%.2f" % (c, percentage),
					horizontalalignment='center',
					verticalalignment='center',
					fontsize=10
				)

			sum_width += percentage


	legend_patches = [
		mpatches.Patch(color=color_dict[c], label=c) for c in classes
	]

	plt.legend(
		handles=legend_patches,
		loc=2,
		bbox_to_anchor=(1.05, 1.05),
		prop={'size':10}
	)

	plt.savefig(img_fname, bbox_inches='tight')

	print("Done")

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--indel_class_fnames", action="append")
	parser.add_argument("-l", "--indel_file_labels", action="append")
	parser.add_argument("plot_title")
	parser.add_argument("img_fname")
	args = parser.parse_args()

	classes = [
		"LINE",
		"SINE",
		"LTR",
		"DNA",
		"Low_complexity",
		"Simple_repeat",
		"Satellite",
		"pseudogene",
		"RNA",
		"Other",
		"Mixed",
		"NONE"
	]

	proportions_dict, totals = get_proportions_dict(
		args.indel_class_fnames,
		args.indel_file_labels,
		classes
	)

	draw_barchart(
		proportions_dict,
		totals,
		args.indel_file_labels,
		classes,
		args.plot_title,
		args.img_fname
	)
