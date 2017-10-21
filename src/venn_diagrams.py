'''
Given 2 or 3 txt files with one item per line, create a Venn diagram of the
items
'''

import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

def draw_venns(input_list, label_list, img_fname):

    if len(input_list) == 2:
        venn_fn = venn2
    elif len(input_list) == 3:
        venn_fn = venn3
    else:
        print("Incorrect number of inputs - must be either 2 or 3")
        sys.exit()

    print(zip(input_list, label_list))

    set_dict = {
        l:set([line.strip() for line in open(fn, 'r')])
        for l,fn in zip(label_list, input_list)}

    fig = plt.figure(figsize=(8,8))
    v = venn_fn([set_dict[l] for l in label_list], set_labels=label_list)
    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", action="append")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    draw_venns(args.input, args.label, args.img_fname)
