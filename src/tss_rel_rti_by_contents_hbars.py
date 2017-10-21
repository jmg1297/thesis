import re
import random as rand
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_plot_data(intersect_fname):
    '''
    Parse the output from bedtools intersect, run on RTI consensus transcripts
    and classified RTIs. Create a dictionary with the following structure:
    RTI classification
        -> total number of RTIs with this classification:<value>
        -> percentage of RTIs with this classification:<value>
        -> number of transcripts with TSS inside RTI:<value>
        -> number of transcripts with TSS outside RTI:<value>
        -> percentage of transcripts with TSS inside:<value>
        -> percentage of transcripts with TSS outside:<value>
    '''
    plot_data = {}
    with open(intersect_fname, 'r') as f:
        for line in f:
            line_list = re.split("\s+", line.strip())
            classification = line_list[12]
            try:
                plot_data[classification]["total"] += 1
            except KeyError:
                plot_data[classification] = {
                    "total":1,
                    "percent":0.0,
                    "num_in":0,
                    "num_out":0,
                    "pc_in":0.0,
                    "pc_out":0.0
                }
            if line_list[5] == "+":
                tss = int(line_list[1])
            elif line_list[5] == "-":
                tss = int(line_list[2])
            else:
                tss = int(line_list[1])

            rti_start = int(line_list[7])
            rti_end = int(line_list[8])

            if rti_start <= tss <= rti_end:
                plot_data[classification]["num_in"] += 1
            else:
                plot_data[classification]["num_out"] += 1

    total_rtis = sum([plot_data[c]["total"] for c in plot_data.keys()])
    for c in plot_data.keys():
        plot_data[c]["percent"] = 100*float(plot_data[c]["total"])/total_rtis
        plot_data[c]["pc_in"] = 100*float(plot_data[c]["num_in"])/plot_data[c]["total"]
        plot_data[c]["pc_out"] = 100*float(plot_data[c]["num_out"])/plot_data[c]["total"]

    return plot_data

def draw_barchart(plot_data, img_fname):
    classifications = sorted(
        plot_data.keys(),
        key=lambda x: plot_data[x]["total"],
        reverse=False
    )

    num_bars = len(classifications)
    height = 0.8

    fig, ax  = plt.subplots()

    plt.suptitle("RTI Transcript TSS Relative Positions")
    ax.set_ylabel("RTI Content Classification")
    ax.set_xlabel("Percentage of Transcript TSSs")

    ax.set_xlim([-5, 105])
    ax.set_ylim([0.5, num_bars+1])
    ax.set_yticks([i+height/2 for i in range(1, num_bars+1)])
    ytick_labels = ["{} ({})".format(c, plot_data[c]["total"]) for c in classifications]
    ax.set_yticklabels(ytick_labels)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    for i,c in enumerate(classifications):
        pc_in = plot_data[c]["pc_in"]
        pc_out = plot_data[c]["pc_out"]

        ax.barh(i+1, pc_in, height=height, color='purple', lw=0, alpha=0.7)
        ax.barh(i+1, pc_out, height=height, color = 'green', lw=0, alpha=0.7, left=pc_in)

    legend_patches = [
        mpatches.Patch(color="purple", label="TSS inside RTI", alpha=0.7),
        mpatches.Patch(color="green", label="TSS outside RTI", alpha=0.7)
    ]

    plt.legend(
        handles=legend_patches, loc=2,
        bbox_to_anchor=(1.05, 1.05), prop={'size':10})

    plt.savefig(img_fname, bbox_inches='tight')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("intersect_fname")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(args.intersect_fname)
    draw_barchart(plot_data, args.img_fname)
