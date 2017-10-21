'''
Script to get the coordinates of indels inside repeats using
rmskJoinedBaseline.txt type file. Extracts all the blocks and their types
(aligned, indel, etc.), then finds gaps just using the aligned blocks, rather
than relying on insertion/deletion starts and sizes
'''

import sys
import re
import subprocess
import argparse
from math import *

def line_to_blocks(line_list):
	'''
	Parse a line from rmskJoinedBaseline.txt to produce a list of the element's
	component blocks and their types
	'''

	# Define a list of dictionaries, one dict per block
	# Keys are start, end, type
	# Options for type are ("aligned", "insertion", "deletion", "indel")
	block_list = []
	num_blocks = int(line_list[10])
	# Get the block sizes and start positions
	block_lengths = [int(x) for x in line_list[11].split(",")]
	block_starts = [int(x) for x in line_list[12].split(",")]

	# The beginning of the whole element, including unaligned sequence at the
	# beginning
	element_start = int(line_list[2])
	element_end = int(line_list[3])

	# The beginning of the aligned sequence
	aligned_start = int(line_list[7])
	aligned_end = int(line_list[8])

	current_coord = element_start

	for i, (size, start) in enumerate(zip(block_lengths, block_starts)):
		if i == 0:
			# First block is always unaligned (maybe empty)
			block_list.append({"start":element_start, "end":aligned_start, "type":"indel"})
			current_coord = aligned_start
		elif i == num_blocks-1:
			# Last block is always unaligned (maybe empty)
			block_list.append({"start":aligned_end, "end":element_end, "type":"indel"})
		else:
			if size > 0 and start >= 0:
				# This is an aligned block
				block_list.append({"start":current_coord, "end":current_coord+size, "type":"aligned"})
				current_coord += size
			elif size > 0 and start == -1:
				# This is an indel or deletion.
				# Size here gives amount of missing sequence from the consensus,
				# rather than the size of the insertion. Need to use the start
				# of the next block to infer the actual size of the insertion

				indel_end = element_start + block_starts[i+1]
				if current_coord == indel_end:
					# Pure deletion with no insertion
					block_list.append({"start":current_coord, "end":indel_end, "type":"deletion"})
				else:
					block_list.append({"start":current_coord, "end":indel_end, "type":"indel"})
					current_coord = indel_end

			elif size <= 0 and start == -1:
				# This is an insertion (without deletion).
				# We can't get the end of the block using the size of the block
				# so we have to use relative start of the *next* block.
				insert_end = element_start + block_starts[i+1]
				block_list.append({"start":current_coord+size, "end":insert_end, "type":"insertion"})
				current_coord = insert_end

	if len(block_list) != num_blocks:
		print("Incorrect number of blocks found")
		print("Expected {}, found {}".format(num_blocks, len(block_list)))
		print("Problem line:")
		print("\t".join(line_list))
		print("Exiting")
		sys.exit()

	return block_list

def get_blocks(rmsk_input):
	'''
	Parse a rmskJoinedBaseline.txt file produce a dict where keys are element
	id numbers and values are lists of component blocks plus other information
	about the element
	'''
	blocks_by_element = {}

	with open(rmsk_input, 'r') as f:

		num_lines = int(subprocess.check_output("wc -l %s | awk '{print $1}'" % rmsk_input, shell = True))
		percent = 0.0

		print("Processing rmsk file")
		print('['+' '*50+']'),
		print('\r'),

		for i,line in enumerate(f):
			line_list = re.split("\s+", line.strip())

			id_num = int(line_list[-1])

			blocks_by_element[id_num] = {
				"chrom":line_list[1],
				"start":int(line_list[2]),
				"end":int(line_list[3]),
				"name":line_list[4],
				"strand":line_list[6],
				"blocks":line_to_blocks(line_list)
			}

			percent = 100.0*(float(i)/num_lines)
			num_bars = int(floor(percent/2))
			print('['+'|'*num_bars+' '*(50-num_bars)+']'),
			print('\r'),
			sys.stdout.flush()

		print('['+'|'*50+']')
		print('\r'),
		print("Done")
		sys.stdout.flush()

	return blocks_by_element

def find_gaps(block_list):
	'''
	Function to find the gaps using the positions of the aligned blocks
	'''
	current_start = -1
	current_end = -1
	gap_list = []
	for block in block_list:
		if block["type"] == "aligned":
			if current_start == -1 and current_end == -1:
				current_start = block["start"]
				current_end = block["end"]
			else:
				if block["start"]-current_end > 1:
					gap_list.append({"start":current_end, "end":block["start"]})
				current_start = block["start"]
				current_end = block["end"]
	return gap_list

def write_gaps(blocks_by_element, output_fname):
	'''
	Extract the indel blocks and write them to the given file
	'''

	print("Extracting insertions and indels")

	output_list = []
	num_lost = 0

	print('['+' '*50+']'),
	print('\r'),
	num_elems = len(blocks_by_element)
	percent = 0.0

	for i,id_num in enumerate(blocks_by_element.keys()):
		element_info = blocks_by_element[id_num]
		name = element_info["name"]+":"+str(id_num)

		gap_list = find_gaps(element_info["blocks"])
		for gap in gap_list:
			out_str = "\t".join([
				element_info["chrom"],
				str(gap["start"]),
				str(gap["end"]),
				name,
				"0",
				element_info["strand"]
			])
			output_list.append(out_str)

		percent = 100.0*(float(i)/num_elems)
		num_bars = int(floor(percent/2))
		print('['+'|'*num_bars+' '*(50-num_bars)+']'),
		print('\r'),
		sys.stdout.flush()

	print('['+'|'*50+']')
	print("Done")
	print("Writing to file ... "),
	sys.stdout.flush()
	with open(output_fname, 'wa') as f:
		f.write("\n".join(output_list))
	print("Done")


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("rmsk_input")
	parser.add_argument("output_fname")
	args = parser.parse_args()

	blocks_by_element = get_blocks(args.rmsk_input)
	write_gaps(blocks_by_element, args.output_fname)































#
