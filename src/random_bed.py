'''
Script to randomly choose regions from one file (the query file) according to a length
distribution defined by another set of regions (the reference file).
'''

import sys
import re
import numpy.random as rand
import linecache
import subprocess
from math import *

def make_length_dict(query_bed_fname):
	'''
	Go through the query file and build up a dictionary with region lengths as
	keys and lists of BED-format region lines as values
	'''

	print("Creating length dictionary ..."),
	length_dict = {}
	with open(query_bed_fname, 'r') as f:
		for line in f:
			line_list = re.split("\s+", line.strip())
			length = int(line_list[2]) - int(line_list[1])
			try:
				length_dict[length].append(line)
			except KeyError:
				length_dict[length] = [line]
	print("Done")
	return length_dict

def pick_random(length_dict, ref_length):
	'''
	Get a random region from the query file based on a length from the reference
	file
	'''
	try:
		# If there are regions with exactly the right length, take this list
		region_list = length_dict[ref_length]
	except KeyError:
		# If there are no regions of exactly the right length, find the list of regions with
		# the closest length
		distances = [[abs(ref_length - key_length), key_length] for key_length in sorted(length_dict.keys())]
		ref_length = min(distances)[1]
		region_list = length_dict[ref_length]

	# Choose a region at random from the list and then remove it (sampling without replacement)
	rand_region = rand.choice(region_list)
	region_list.remove(rand_region)

	if len(region_list) > 0:
		length_dict[ref_length] = region_list
	else:
		# If removing the region left an empty list, delete that key-value pair from the dictionary
		removed = length_dict.pop(ref_length, None)

	# Return the random region and the changed lengthDict
	return rand_region, length_dict

def get_sample(ref_bed_fname, query_bed_fname, out_bed_fname, jitter_length, sample_size):
	# Get the number of lines in the reference file
	num_ref_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'" % ref_bed_fname, shell = True))

	# Make the length dictionary
	length_dict = make_length_dict(query_bed_fname)

	out_lines = []
	for i in range(sample_size):
		ref_line_num = rand.randint(1, num_ref_lines+1)
		ref_line = linecache.getline(ref_bed_fname, ref_line_num)
		ref_line_list = re.split("\s+", ref_line)
		ref_length = int(ref_line_list[2]) - int(ref_line_list[1])

		# To prevent spikes at exactly the same lengths as in the reference file, adjust the
		# lengths by a small amount
		jittered_ref_length = rand.randint(ref_length - jitter_length, ref_length + jitter_length)

		# Get a random region
		rand_region, length_dict = pick_random(length_dict, jittered_ref_length)
		out_lines.append(rand_region)

	# Write the random regions to the output file
	with open(out_bed_fname, 'w') as f:
		f.write("".join(out_lines))

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("ref_bed_fname")
	parser.add_argument("query_bed_fname")
	parser.add_argument("out_bed_fname")
	parser.add_argument("jitter_length", type=int)
	parser.add_argument("sample_size", type=int)
	args = parser.parse_args()

	get_sample(
		args.ref_bed_fname,
		args.query_bed_fname,
		args.out_bed_fname,
		args.jitter_length,
		args.sample_size
	)
