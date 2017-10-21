'''
Script to produce volcano plots of fold changes for three sets of transcripts.
One set is B cell-specific, one set is T cell-specific, and one is shared by
both cell types.
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def read_ballgown(ballgown_results_fname):
    ballgown_dict = {}
    with open(ballgown_results_fname, 'r') as f:
        keys = f.readline()
        for line in f:
            if "NA" in line:
                continue
            line_list = line.strip().split()
            transcript_id = line_list[0]
            fc, pval, qval = map(float, line_list[1:])
            ballgown_dict[transcript_id] = {"fc":fc, "pval":pval, "qval":qval}
    return ballgown_dict

def get_plot_values(transcript_list_fname, ballgown_dict):
    print(transcript_list_fname)

    with open(transcript_list_fname, 'r') as f:
        ref_transcripts = [l.strip() for l in f]

    plot_values = []
    for t in ref_transcripts:
        try:
            x = np.log2(ballgown_dict[t]["fc"])
            y = -np.log10(ballgown_dict[t]["pval"])
            if ballgown_dict[t]["pval"] < 0.05:
                print([t, x])
        except KeyError:
            continue
        else:
            plot_values.append([x,y])

    return plot_values

def get_plot_data(B_spec_parents, T_spec_parents, shared_parents, ballgown_results_fname):
    ballgown_dict = read_ballgown(ballgown_results_fname)

    plot_data = {
        "B_spec":get_plot_values(B_spec_parents, ballgown_dict),
        "T_spec":get_plot_values(T_spec_parents, ballgown_dict),
        "shared":get_plot_values(shared_parents, ballgown_dict)
    }

    return plot_data

def draw_volcano(plot_data, img_fname, plot_title):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,8)

    color_dict = {"B_spec":"blue", "T_spec":"red", "shared":"grey"}
    size_dict = {"B_spec":20, "T_spec":20, "shared":6}
    alpha_dict = {"B_spec":0.7, "T_spec":0.7, "shared":0.4}
    label_dict = {"B_spec":"B-specific", "T_spec":"T-specific", "shared":"Both\ncell types"}

    for k in ["shared", "T_spec", "B_spec"]:
        data = plot_data[k]
        x = []
        y = []
        for t in data:
            x.append(t[0])
            y.append(t[1])
        ax.scatter(
            x, y,
            lw=0, color=color_dict[k], s=size_dict[k], alpha=alpha_dict[k],
            label=label_dict[k]
        )

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

    plt.legend(bbox_to_anchor=(0., 0.2), loc=0, borderaxespad=0., fontsize=10)

    plt.savefig(img_fname)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("B_spec")
    parser.add_argument("T_spec")
    parser.add_argument("shared")
    parser.add_argument("ballgown_results_fname")
    parser.add_argument("img_fname")
    parser.add_argument("plot_title")
    args = parser.parse_args()

    plot_data = get_plot_data(
        args.B_spec, args.T_spec, args.shared,
        args.ballgown_results_fname
    )

    draw_volcano(plot_data, args.img_fname, args.plot_title)
