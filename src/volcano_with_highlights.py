'''
Draw a volcano plot from a Ballgown-derived TSV file, highlighting genes from a
given list
'''

import matplotlib.pyplot as plt
import numpy as np

def read_ballgown(ballgown_fname):
    ballgown_dict = {}
    with open(ballgown_fname, 'r') as f:
        keys = f.readline()
        for line in f:
            if "NA" in line:
                continue
            line_list = line.strip().split()
            transcript_id = line_list[0]
            fc, pval, qval = map(float, line_list[1:])
            ballgown_dict[transcript_id] = {"fc":fc, "pval":pval, "qval":qval}
    return ballgown_dict

def read_transcript_list(transcript_list_fname):
    with open(transcript_list_fname, 'r') as f:
        transcripts = [l.strip() for l in f]
    return transcripts

def draw_volcano(ballgown_dict, highlights, img_fname, plot_title):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,8)

    hlx = []
    hly = []

    bgx = []
    bgy = []

    for transcript, data in ballgown_dict.iteritems():
        if transcript in highlights:
            hlx.append(np.log2(data["fc"]))
            hly.append(-np.log10(data["pval"]))
        else:
            bgx.append(np.log2(data["fc"]))
            bgy.append(-np.log10(data["pval"]))

    ax.scatter(hlx, hly, lw=0, c="red", s=20, alpha=0.7, zorder=1)
    ax.scatter(bgx, bgy, lw=0, c="gray", s=6, alpha=0.4, zorder=0)

    ax.axhline(y=-np.log10(0.05), color="black", linestyle="dashed")
    ax.text(-6, -np.log10(0.05)*1.15, s="$p=0.05$", color="black")

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.set_xlim([-abs(max(map(abs, ax.get_xlim()))), abs(max(map(abs, ax.get_xlim())))])
    #ax.set_ylim([0,12])

    ax.set_xlabel("$\log_2(FC)$")
    ax.set_ylabel("$-\log_{10}(p)$", rotation="horizontal")
    ax.yaxis.set_label_coords(0.5, 1.05)

    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("ballgown_fname")
    parser.add_argument("transcript_list_fname")
    parser.add_argument("img_fname")
    parser.add_argument("plot_title")
    args = parser.parse_args()

    ballgown_dict = read_ballgown(args.ballgown_fname)
    highlights = read_transcript_list(args.transcript_list_fname)
    draw_volcano(ballgown_dict, highlights, args.img_fname, args.plot_title)
