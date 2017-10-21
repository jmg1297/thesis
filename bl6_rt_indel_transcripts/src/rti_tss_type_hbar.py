'''
Given labelled mBED files with transcripts and corresponding RTIs, for each
file, for each transcript, decide whether the TSS is inside the RTI, inside the
upstream host block, or upstream of the upstream host block. Plot results in a
horizontal barchart.
'''

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_rt_block_dict(rt_block_fname):
    rt_block_dict = {}
    with open(rt_block_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            rt_num = line_list[3].split(":")[3]
            start = int(line_list[1])
            end = int(line_list[2])
            try:
                rt_block_dict[rt_num].append({"start":start, "end":end})
            except KeyError:
                rt_block_dict[rt_num] = [{"start":start, "end":end}]
    return rt_block_dict

def get_upstream_block(rti_start, rti_end, rti_name, strand, rt_block_dict):
    rti_num = rti_name.split(":")[3]
    all_blocks = rt_block_dict[rti_num]
    for block in all_blocks:
        if strand == "+":
            if rti_start == block["end"]:
                return block["start"], block["end"]
        else:
            if rti_end == block["start"]:
                return block["start"], block["end"]

def process_mbed(mbed_fname, rt_block_dict):
    counts = {"Outside Host":0, "Inside Host":0, "Inside RTI":0}

    with open(mbed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            strand = line_list[5]
            if strand == "+":
                tss = int(line_list[1])
            else:
                tss = int(line_list[2])

            rti_start = int(line_list[7])
            rti_end = int(line_list[8])

            if rti_start < tss < rti_end:
                counts["Inside RTI"] += 1
            else:
                rti_name = line_list[9]
                host_start, host_end = get_upstream_block(rti_start, rti_end, rti_name, strand, rt_block_dict)
                if host_start <= tss <= host_end:
                    counts["Inside Host"] += 1
                else:
                    counts["Outside Host"] += 1
    total = sum(counts.values())

    props = {k:100.0*float(v)/total for k,v in counts.items()}

    return props, total

def get_plot_data(mbed_list, label_list, rt_block_fname):
    rt_block_dict = get_rt_block_dict(rt_block_fname)

    plot_data = {}
    totals = {}

    for fn,l in zip(mbed_list, label_list):
        plot_data[l], totals[l] = process_mbed(fn, rt_block_dict)

    return plot_data, totals

def draw_hbars(plot_data, totals, label_list, img_fname):
    color_dict = {"Outside Host":"red", "Inside Host":"blue", "Inside RTI":"green"}
    classes = ["Outside Host", "Inside Host", "Inside RTI"]
    num_bars = len(label_list)

    fig = plt.figure(figsize = [12,3*num_bars])

    axis = fig.add_axes([0.05,0.05,0.9,0.9])

    axis.set_ylabel("Region set")
    axis.set_xlabel("Percentage of regions")

    height = 0.8

    axis.set_xlim([-5,105])
    axis.set_ylim([0.5, num_bars+1])
    axis.set_yticks([i + height/2 for i in range(1, num_bars+1)])
    ytick_labels = [
        "{} ({})".format(label, totals[label]) for label in label_list]
    axis.set_yticklabels(ytick_labels)

    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    for i,label in enumerate(label_list):
        data = plot_data[label]
        sum_width = 0
        for c in classes:
            percentage = data[c]
            axis.barh(
                i+1,
                percentage,
                height=height,
                left=sum_width,
                color=color_dict[c],
                lw=0,
                alpha=0.7
            )

            axis.text(
                sum_width + percentage/2,
                i + 1 + height/2,
                "%s\n%.2f" % (c, percentage),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=10
            )

            sum_width += percentage

    legend_patches = [
        mpatches.Patch(color=color_dict[c], label=c) for c in classes]

    plt.legend(
        handles=legend_patches,
        loc=2,
        bbox_to_anchor=(1.05, 1.05),
        prop={'size':10}
    )

    plt.savefig(img_fname, bbox_inches='tight')

    print("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mbed", action="append")
    parser.add_argument("-l", "--label", action="append")
    parser.add_argument("rt_block_fname")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    plot_data, totals = get_plot_data(args.mbed, args.label, args.rt_block_fname)

    draw_hbars(plot_data, totals, args.label, args.img_fname)
