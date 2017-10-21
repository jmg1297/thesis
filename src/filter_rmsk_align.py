import sys
import subprocess
from StringIO import StringIO
import re

intersect_fname = sys.argv[1]

allowed = [
    line.strip() \
        for line in StringIO(
            subprocess.check_output(
                "cat %s | awk '{print $10}' | cut -d: -f1,2,3 | sort | uniq | sed '1d'" % intersect_fname,
                shell=True
            )
        ).readlines()
]

rmsk_align_bed_fname = sys.argv[2]
output_fname = sys.argv[3]

with open(rmsk_align_bed_fname, 'r') as f, open(output_fname, 'wa') as out:
    for line in f:
        element = ":".join(re.split("\s+", line.strip())[3].split(":")[0:-1])
        if element in allowed:
            out.write(line)
