'''
Use a GTF output by StringTie to create a database with two tables, transcripts
and exons.
'''

import sqlite3
import re
from collections import OrderedDict

def scrub(s):
    return re.sub("[\";]", "", s)

def gtf_to_db(gtf_fname):
    db_fname = gtf_fname[0:-4] + ".db"

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    cur.execute(
        '''
        CREATE TABLE transcript (
            chromosome text,
            start_coord integer,
            end_coord integer,
            strand text,
            gene_id text,
            transcript_id text,
            num_exons integer,
            exon_length integer,
            PRIMARY KEY (transcript_id)
        )
        '''
    )

    cur.execute(
        '''
        CREATE TABLE exon (
            chromosome text,
            start_coord integer,
            end_coord integer,
            strand text,
            gene_id text,
            transcript_id text,
            exon_number integer,
            PRIMARY KEY (transcript_id, exon_number)
        )
        '''
    )

    transcript_dict = OrderedDict()
    exon_data = []

    with open(gtf_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            line_list = line.strip().split()
            table = line_list[2]
            chromosome = line_list[0]
            if "chr" not in chromosome:
                chromosome = "chr" + chromosome
            start_coord = int(line_list[3])
            end_coord = int(line_list[4])
            strand = line_list[6]
            gene_id = scrub(line_list[9])
            transcript_id = scrub(line_list[11])
            if table == "exon":
                exon_number = int(scrub(line_list[13]))
                exon_data.append(
                    (
                        chromosome, start_coord, end_coord, strand,
                        gene_id, transcript_id, exon_number
                    )
                )
                transcript_dict[transcript_id]["exon_length"] \
                    += (end_coord - start_coord)
                transcript_dict[transcript_id]["num_exons"] += 1
            else:
                transcript_dict[transcript_id] = OrderedDict([
                    ("chromosome", chromosome),
                    ("start_coord", start_coord),
                    ("end_coord", end_coord),
                    ("strand", strand),
                    ("gene_id", gene_id),
                    ("transcript_id", transcript_id),
                    ("num_exons", 0),
                    ("exon_length", 0)
                ])

    transcript_data = [d.values() for d in transcript_dict.values()]

    tcols = len(transcript_data[0])
    ecols = len(exon_data[1])

    cur.executemany(
        "INSERT INTO transcript VALUES ({})".format(",".join(["?"]*tcols)),
        transcript_data
    )

    cur.executemany(
        "INSERT INTO exon VALUES ({})".format(",".join(["?"]*ecols)),
        exon_data
    )

    conn.commit()

    num_transcripts = cur.execute("SELECT COUNT(*) FROM transcript").fetchone()
    num_exons = cur.execute("SELECT COUNT(*) FROM exon").fetchone()

    print("{} transcripts and {} exons inserted into database {}"
            .format(num_transcripts[0], num_exons[0], db_fname))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_fname")
    args = parser.parse_args()

    gtf_to_db(args.gtf_fname)
