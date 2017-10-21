'''
Create a JSON file of BL6 sample GTF filenames and colors
'''

import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_gtf_json(sample_dir, gtf_fname, output_fname):
    gtf_dict = {}
    cmap = plt.get_cmap('tab20c')
    colors = [cmap(i/20.0) for i in range(16) if (i+1)%4 != 0]
    color_list = [cmap(i/20.0) for i in range(16) if (i+1)%4 != 0]
    for i,sample in enumerate(sorted(os.listdir(sample_dir))):
        fpath = os.path.join(sample_dir, sample, gtf_fname)
        color = color_list[i]
        gtf_dict[sample] = {"fname":fpath, "color":color}
    with open(output_fname, 'wa') as f:
        json.dump(gtf_dict, f)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_dir")
    parser.add_argument("fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    write_gtf_json(args.sample_dir, args.fname, args.output_fname)
