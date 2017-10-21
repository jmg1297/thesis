'''
Script to classify retrotransposon indels based on their contents.
Input is the output from bedtools intersect. Output is a BED-like file with two
extra columns.
First extra column is the classification. This is the intersection with the
greatest percentage overlap, providing that the percentage is above a given
threshold.
Second extra column contains information about all the intersections, in the
following format:
<class>:<percentage>,<class>:<percentage>,...
e.g., LINE:80,Simple_repeat:15
'''

import sys
import re
import subprocess
from math import *
import argparse

class Indel(object):
	def __init__(self, start, end, strand, min_percent_overlap):
		self.length = int(end) - int(start)
		self.strand = strand
		self.min_percent_overlap = float(min_percent_overlap)
		self.intersections = []

	def add_intersection(self, intersect_str, overlap):
		if overlap != 0:
			overlap_percent = 100*(float(overlap)/self.length)
			self.intersections.append([intersect_str, overlap_percent])
		else:
			self.intersections.append(["NONE", 100.0])

	def get_intersect_list(self):
		intersect_dict = {}
		for intersect_str, overlap_percent in self.intersections:
			if intersect_str != "NONE":
				elem_info = intersect_str.split(" ")
				type_info = elem_info[3]
				elem_type = type_info.split(":")[0]
			else:
				elem_type = "NONE"

			try:
				intersect_dict[elem_type] += overlap_percent
			except KeyError:
				intersect_dict[elem_type] = overlap_percent

		total_percent = sum(intersect_dict.values())
		if total_percent < 100:
			intersect_dict["NONE"] = 100 - total_percent

		intersect_list = sorted(
			[
				[intersect_dict[elem_type], elem_type] \
				for elem_type in intersect_dict.keys()
			],
			reverse=True
		)

		major_intersect = intersect_list[0]
		if major_intersect[0] >= self.min_percent_overlap:
			classification = major_intersect[1]
		else:
			classification = "Mixed"

		return intersect_list, classification

def intersect_file_to_dict(intersect_fname, min_percent_overlap):
	'''
	Keys are a string with chromosome, start, end, name for an indel.
	Values are instances of the Indel class.
	'''

	print("Processing intersection file ... ")
	sys.stdout.flush()

	intersect_dict = {}

	with open(intersect_fname, 'r') as f:

		print('['+' '*50+']'),
		print('\r'),
		num_lines = int(
			subprocess.check_output(
				"wc -l %s | awk '{print $1}'" % intersect_fname,
				shell=True
			)
		)
		percent = 0.0

		for i,line in enumerate(f):
			line_list = re.split("\s+", line.strip())
			indel_name = " ".join(line_list[0:4])
			intersect_str = " ".join(line_list[6:10])
			try:
				intersect_dict[indel_name].add_intersection(
					intersect_str,
					int(line_list[-1])
				)
			except KeyError:
				intersect_dict[indel_name] = Indel(
											line_list[1],
											line_list[2],
											line_list[5],
											min_percent_overlap
				)
				intersect_dict[indel_name].add_intersection(
					intersect_str,
					int(line_list[-1])
				)

			percent = 100.0*(float(i)/num_lines)
			num_bars = int(floor(percent/2))
			print('\r'),
			print('['+'|'*num_bars+' '*(50-num_bars)+']'),
			sys.stdout.flush()

	print("\nDone")
	return intersect_dict

def write_output(intersect_dict, output_fname):
	print("Writing output file ... "),
	sys.stdout.flush()

	with open(output_fname, 'wa') as f:
		for indel_name in intersect_dict.keys():
			out_list = [
				re.sub(" ", "\t", indel_name),
				"0",
				intersect_dict[indel_name].strand
			]

			intersect_list, classification \
				= intersect_dict[indel_name].get_intersect_list()

			intersect_str = ",".join(
				[
					"{}:{}".format(elem_type, percent) \
					for percent, elem_type in intersect_list
				]
			)
			out_list += [classification, intersect_str]
			out_str = "\t".join(out_list) + "\n"
			f.write(out_str)

	print("Done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("intersect_fname")
	parser.add_argument("output_fname")
	parser.add_argument("min_percent_overlap")
	args = parser.parse_args()

	intersect_dict = intersect_file_to_dict(
		args.intersect_fname,
		args.min_percent_overlap
	)
	write_output(intersect_dict, args.output_fname)

















#
