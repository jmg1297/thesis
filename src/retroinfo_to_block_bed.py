# Script to take a ucscRetroInfo6.txt formatted file and output a BED file
# containing the blocks ("exons") of each retrogene

import sys
import re
import subprocess
from math import *
import argparse

def get_blocks(line, line_num):
	'''
	Extract a list of blocks from a retrogene line, then return the lines to
	write to the output file with formatted names
	'''

	line_list = re.split("\s+", line.strip())
	chromosome = line_list[0]
	start = int(line_list[1])
	end = int(line_list[2])
	name = line_list[3]
	strand = line_list[5]
	block_sizes = [int(x) for x in line_list[10].split(",") if x != '']
	block_rel_starts = [int(x) for x in line_list[11].split(",") if x != '']
	output_lines = []
	for i, (st, sz) in enumerate(zip(block_rel_starts, block_sizes)):
		block_start = start + st
		block_end = block_start + sz
		block_name = ":".join([
			"pseudogene",
			"pseudogene",
			name,
			str(line_num)+"."+str(i+1)
		])

		output_lines.append(
			"\t".join([
				chromosome,
				str(block_start),
				str(block_end),
				block_name,
				"0",
				strand
			])
		)
	return "\n".join(output_lines) + "\n"

def retroinfo_to_bed(retroinfo_fname, out_fname):
	print("Beginning block extraction ... ")
	with open(retroinfo_fname, 'r') as f, open(out_fname, 'wa') as out:
		print('['+' '*50+']'),
		print('\r'),
		num_lines = int(
			subprocess.check_output(
				"wc -l %s | awk '{print $1}'" % retroinfo_fname,
				shell =True
			)
		)
		percent = 0.0
		i = 1

		for line in f:
			out.write(get_blocks(line, i))
			i += 1
			percent = 100.0*(float(i-1)/num_lines)
			num_bars = int(floor(percent/2))
			print('\r'),
			print('['+'|'*num_bars+' '*(50-num_bars)+']'),
			sys.stdout.flush()

	print("\nDone")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("retroinfo_fname")
	parser.add_argument("out_fname")
	args = parser.parse_args()

	retroinfo_to_bed(args.retroinfo_fname, args.out_fname)
