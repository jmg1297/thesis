'''
Create a DB of Ensembl transcripts from a GTF file
'''

import sys
import sqlite3
import re

def scrub(s):
    return re.sub("[\";]", "", s)

def get_attr_dict(attr_list):
    N = len(attr_list)

    keys = [attr_list[i] for i in range(0, N, 2)]
    vals = [scrub(attr_list[i]) for i in range(1, N, 2)]

    return {k:v for k,v in zip(keys, vals)}

def make_transcript_table(input_fname, db_fname):

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS transcript")

    cur.execute(
        '''
        CREATE TABLE transcript (
        chromosome text,
        start_coord integer,
        end_coord integer,
        strand text,
        gene_id text,
        transcript_id text,
        gene_name text,
        gene_biotype text,
        transcript_biotype text
        )
        '''
    )

    data = []

    with open(input_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split()

            if line_list[2] != "transcript":
                continue

            chromosome = line_list[0]
            start_coord = int(line_list[3])
            end_coord = int(line_list[4])
            strand = line_list[6]
            attr_dict = get_attr_dict(line_list[8:])

            data.append((
                chromosome, start_coord, end_coord, strand,
                attr_dict["gene_id"], attr_dict["transcript_id"],
                attr_dict["gene_name"],
                attr_dict["gene_biotype"], attr_dict["transcript_biotype"]
            ))

    num_cols = len(data[0])
    print(num_cols)
    cur.executemany(
        "INSERT INTO transcript VALUES ({})".format(",".join(['?']*num_cols)),
        data
    )

    conn.commit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fname")
    parser.add_argument("db_fname")
    args = parser.parse_args()

    make_transcript_table(args.input_fname, args.db_fname)
