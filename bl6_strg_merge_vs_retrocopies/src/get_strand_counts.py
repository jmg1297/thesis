'''
Create a JSON file of counts of different strand combinations
'''

import sys
import subprocess
import json

input_mbed = sys.argv[1]
output_json = sys.argv[2]

strands = ["+", "-"]
counts_dict = {s:{t:{u:0 for u in strands} for t in strands} for s in strands}

cmd = "cat %s | awk '{print $6,$12,$19}' | sort | uniq -c | sed '/\./d'" % input_mbed
raw_output = subprocess.check_output(cmd, shell=True)

for line in raw_output.strip().split("\n"):
    count,  transcript_strand, rc_strand, parent_strand = line.split()
    counts_dict[parent_strand][rc_strand][transcript_strand] = int(count)

with open(output_json, 'wa') as out:
    json.dump(counts_dict, out)
