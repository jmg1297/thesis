'''
Given a set of mBED files with transcripts, retrocopies, and parent transcripts,
find the retrocopies that are specific to each. Write results to files.
'''

import tempfile as tf
import os
import subprocess

def get_rc_dict(filename):
    rc_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            retrocopy = line.strip().split()[9]
            try:
                rc_dict[retrocopy].append(line)
            except KeyError:
                rc_dict[retrocopy] = [line]

    return rc_dict

def get_spec_retrocopies(input_fname_list):

    input_dicts = {fn:get_rc_dict(fn) for fn in input_fname_list}

    for fn in input_fname_list:
        others = [input_dicts[k].keys() for k in input_fname_list if k != fn]
        spec_retrocopies = set(input_dicts[fn].keys()).difference(*map(set, others))

        out_fn = os.path.join(
            os.path.dirname(fn),
            "spec_expr_retrocopies.mbed"
        )
        out_lines = ""
        for rc in list(spec_retrocopies):
            for line in input_dicts[fn][rc]:
                out_lines += line

        with tf.NamedTemporaryFile() as tmp:
            tmp.write(out_lines)
            tmp.flush()
            subprocess.check_output("sort -k1,1 -k2,2n %s > %s" % (tmp.name, out_fn), shell=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fname", action="append")
    args = parser.parse_args()

    get_spec_retrocopies(args.input_fname)
