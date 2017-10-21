'''
Given two BED files and corresponding files of chromosome sizes, for each
region in the BED file calculate the relative position in the chromosome.
Draw a scatter plot of the relative positions.
'''

import matplotlib.pyplot as plt

def get_chrom_sizes(chrom_size_fname):
    with open(chrom_size_fname, 'r') as f:
        chrom_sizes = {line.strip().split()[0]:int(line.strip().split()[1])
                        for line in f}
    return chrom_sizes

def get_plot_data(ref_label, ref_bed, ref_chrom_size_fname,
        query_label, query_hits, query_chrom_size_fname):

    plot_data = {ref_label:[], query_label:[], "lines":[]}

    ref_chrom_sizes = get_chrom_sizes(ref_chrom_size_fname)
    query_chrom_sizes = get_chrom_sizes(query_chrom_size_fname)

    ref_data = {}

    with open(ref_bed, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            chrom = line_list[0]
            midpoint = (int(line_list[1]) + int(line_list[2]))/2.0
            rel_pos = midpoint/ref_chrom_sizes[chrom]
            name = line_list[3]
            ref_data[name] = rel_pos

    with open(query_hits, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            name = line_list[0]
            chrom = line_list[1]
            if "chr" not in chrom:
                continue
            midpoint = (int(line_list[8]) + int(line_list[9]))/2.0
            rel_pos = midpoint/query_chrom_sizes[chrom]
            plot_data[query_label].append(rel_pos)
            plot_data[ref_label].append(ref_data[name])
            plot_data["lines"].append(line)

    return plot_data

def find_outliers(X, Y, lines):
    upper_c = 0.01
    lower_c = -upper_c
    step = 0.001
    threshold = 0.8
    total = len(X)

    while True:
        count = 0
        close_hits = []
        for x,y,l in zip(X, Y, lines):
            if y + x + lower_c and y < x + upper_c:
                count += 1
                close_hits.append(l)
        if float(count)/total > threshold:
            break
        else:
            upper_c += step
            lower_c = -upper_c

    return upper_c, lower_c, close_hits

def scatter_data(plot_data, ref_label, query_label, plot_title, img_fname, out_fname):

    fig, ax = plt.subplots(1,1)
    ax.scatter(plot_data[ref_label], plot_data[query_label], s=2, c='k', alpha=0.3, lw=0)

    upper_c, lower_c, close_hits = find_outliers(
        plot_data[ref_label], plot_data[query_label], plot_data["lines"]
    )

    ax.plot([0,1], [upper_c, 1+upper_c], 'r--')
    ax.plot([0,1], [lower_c, 1+lower_c], 'r--')

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

    ax.set_aspect('equal')

    print(upper_c)

    ax.set_xlabel(ref_label)
    ax.set_ylabel(query_label)

    plt.savefig(img_fname)

    with open(out_fname, 'wa') as out:
        out.write("".join(close_hits))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_label")
    parser.add_argument("ref_bed")
    parser.add_argument("ref_chrom_size")
    parser.add_argument("query_label")
    parser.add_argument("query_hits")
    parser.add_argument("query_chrom_size")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    parser.add_argument("out_fname")
    args = parser.parse_args()

    plot_data = get_plot_data(
        args.ref_label, args.ref_bed, args.ref_chrom_size,
        args.query_label, args.query_hits, args.query_chrom_size
    )

    scatter_data(
        plot_data, args.ref_label, args.query_label,
        args.plot_title, args.img_fname, args.out_fname
    )
