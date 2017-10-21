'''
Get counts of the strand combinations of
(retrocopy transcript, retrocopy annotation, parent transcript)
'''

import sys
import sqlite3

def get_strand_combos(intersect_fname, retrocopy_db, output_fname):

    conn = sqlite3.connect(retrocopy_db)
    cur = conn.cursor()

    get_strand_query = "SELECT retrogenes_parents.retrogene,ensembl_transcripts.strand FROM retrogenes_parents JOIN ensembl_transcripts ON retrogenes_parents.parents=ensembl_transcripts.transcript_id WHERE retrogenes_parents.retrogene IN ({})"

    retrocopies = []

    transcript_strands = {}
    parent_strands = {}
    retrocopy_strands = {}

    with open(intersect_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            retrocopy = line_list[9]
            retrocopies.append(line_list[9])
            transcript_strands[retrocopy] = line_list[5]
            retrocopy_strands[retrocopy] = line_list[11]


    parent_strands = dict(cur.execute(get_strand_query.format(",".join(["?"]*len(retrocopies))), retrocopies).fetchall())

    strand_combos = []

    for rc in retrocopies:
        try:
            strand_combos.append((
                parent_strands[rc],
                retrocopy_strands[rc],
                transcript_strands[rc]
            ))
        except KeyError:
            continue

    counts = {c:0 for c in set(strand_combos)}
    for c in strand_combos:
        counts[c] += 1

    with open(output_fname, 'wa') as out:
        out.write("\t".join(["parent_strand", "retrocopy_strand", "transcript_strand", "count"]) + "\n")
        for combo, count in counts.iteritems():
            out.write("\t".join(list(combo) + [str(count)]) + "\n")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("intersect_fname")
    parser.add_argument("retrocopy_db")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    get_strand_combos(args.intersect_fname, args.retrocopy_db, args.output_fname)
