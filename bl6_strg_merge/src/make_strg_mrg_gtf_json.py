'''
Create a JSON file of StringTie merge GTF filenames and colors
'''

import json
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_gtf_json(sample_dir, output_fname):
    gtf_dict = {}
    cmap = plt.get_cmap('tab10')
    gtf_fnames = [fn for fn in os.listdir(sample_dir) if ".gtf" in fn]
    for i,fn in enumerate(gtf_fnames):
        fpath = os.path.join(sample_dir, fn)
        color = cmap(i/10.0)
        label = re.split("merge", fn)[0][0:-1]
        gtf_dict[label] = {"fname":fpath, "color":color}
    with open(output_fname, 'wa') as f:
        json.dump(gtf_dict, f, sort_keys=True, indent=4)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_dir")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    write_gtf_json(args.sample_dir, args.output_fname)
