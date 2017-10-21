'''
Given a multibed file of transcripts with corresponding retrocopies, find the
parent transcript for each retrocopy and append the relevant information to the
file in extra columns.
'''

import sqlite3
import subprocess

def make_parent_line(l):
    l.insert(4, 0)
    return "\t".join(map(str, l)) + "\n"

def append_parent_bed(input_fname, rc_name_field, retrocopy_db, output_fname):
    conn = sqlite3.connect(retrocopy_db)
    cur = conn.cursor()

    cmd_str = "cat %s | awk '{print $%s}' | sort | uniq" % (input_fname, rc_name_field)
    retrocopies = subprocess.check_output(cmd_str, shell=True).strip().splitlines()

    query = "SELECT rp.retrogene, et.chromosome, et.start_coord, et.end_coord, et.transcript_id, et.strand \
            FROM retrogenes_parents AS rp \
            JOIN ensembl_transcripts AS et \
            ON rp.parents = et.transcript_id \
            WHERE rp.retrogene IN ({})"\
            .format(",".join(["?"]*len(retrocopies)))

    rows = cur.execute(query, retrocopies).fetchall()

    retro_to_parent = {r[0]:make_parent_line(list(r[1:])) for r in rows}

    with open(input_fname, 'r') as f, open(output_fname, 'wa') as out:
        for line in f:
            retrocopy = line.strip().split()[int(rc_name_field)-1]
            try:
                parent_line = retro_to_parent[retrocopy]
            except KeyError:
                parent_line = "\t".join(["."]*6) + "\n"
            out_line = line.strip() + "\t" + parent_line
            out.write(out_line)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fname")
    parser.add_argument("rc_name_field")
    parser.add_argument("retrocopy_db")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    append_parent_bed(args.input_fname, args.rc_name_field, args.retrocopy_db, args.output_fname)
