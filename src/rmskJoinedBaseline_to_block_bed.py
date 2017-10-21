# Script to go from a rmskJoinedBaseline.txt file to a BED file of aligned blocks

import sys
import re
import subprocess
from math import *
import argparse

def format_name(name_str):
	if "/" in name_str:
		name_list = re.split("[#/:]", name_str)
		return ":".join([
			name_list[1],
			name_list[2],
			name_list[0],
			name_list[3]
		])
	else:
		name_list = re.split("[#:]", name_str)
		return ":".join([
			name_list[1],
			name_list[1],
			name_list[0],
			name_list[2]
		])

def process_line(line):
	'''
	Take a line from the input file and extract the information for a BED file
	'''

	line_list = re.split("\s+", line.strip())
	chrom = line_list[1]
	strand = line_list[6]
	name = format_name(line_list[4]+":"+line_list[-1])
	rel_starts = [int(x) for x in line_list[12].split(",")]
	sizes = [int(x) for x in line_list[11].split(",")]
	abs_start = int(line_list[2])

	out_lines = []

	for start, size, i in zip(rel_starts, sizes, range(len(sizes))):
		if start == -1:
			continue
		else:
			out_lines.append(
				"\t".join(
					[
						chrom,
						str(abs_start + start),
						str(abs_start + start + size),
						name,
						"0",
						strand
					]
				)
			)

	return out_lines

def rmsk_to_bed(rmsk_fname, out_fname):
	'''
	Function to take the name of a rmskJoinedBaseline table file from the UCSC
	Genome Browser and use it to create a BED file of the individual blocks
	'''

	with open(rmsk_fname, 'r') as f:

		print("Processing table file ... ")
		print('['+' '*50+']'),
		print('\r'),
		num_lines = int(
			subprocess.check_output(
				"wc -l %s | awk '{print $1}'" % rmsk_fname,
				shell=True
			)
		)
		percent = 0.0
		i = 0

		out_lines = []

		for line in f:
			out_lines += process_line(line)

			i += 1
			percent = 100.0*(float(i)/num_lines)
			num_bars = int(floor(percent/2))
			print('\r'),
			print('['+'|'*num_bars+' '*(50-num_bars)+']'),
			sys.stdout.flush()

		print("\nDone")

	out_str = "\n".join(out_lines) + "\n"
	with open(out_fname, 'w') as f:
		print("Writing results ... "),
		f.write(out_str)
		print("Done")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("rmsk_fname")
	parser.add_argument("out_fname")
	args = parser.parse_args()

	rmsk_to_bed(args.rmsk_fname, args.out_fname)
