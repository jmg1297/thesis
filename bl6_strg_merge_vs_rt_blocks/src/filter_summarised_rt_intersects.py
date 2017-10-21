'''
Given a BED-like file of transcripts with summarised retrotransposon overlaps,
filter them so that only those with >x% exonic retrotransposon content are
retained
'''

def get_exon_length_dict(exon_bed_fname):
    exon_length_dict = {}
    with open(exon_bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = ".".join(line_list[3].split(".")[0:3])
            length = int(line_list[2]) - int(line_list[1])
            try:
                exon_length_dict[transcript] += length
            except KeyError:
                exon_length_dict[transcript] = length
    return exon_length_dict

def filter_overlaps(summary_bed_fname, exon_bed_fname, overlap_cutoff, output_fname):

    exon_length_dict = get_exon_length_dict(exon_bed_fname)

    with open(summary_bed_fname, 'r') as f, open(output_fname, 'wa') as out:
        for line in f:
            line_list = line.strip().split()
            total_exon_length = exon_length_dict[line_list[3]]
            overlap_total = reduce(
                lambda a,x: a+x,
                map(lambda y: int(y.split(",")[1]), line_list[-1].split(";")),
                0
            )
            if float(overlap_total)/total_exon_length >= overlap_cutoff:
                out.write(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("summary_bed_fname")
    parser.add_argument("exon_bed_fname")
    parser.add_argument("overlap_cutoff", type=float)
    parser.add_argument("output_fname")
    args = parser.parse_args()

    filter_overlaps(args.summary_bed_fname, args.exon_bed_fname, args.overlap_cutoff, args.output_fname)
