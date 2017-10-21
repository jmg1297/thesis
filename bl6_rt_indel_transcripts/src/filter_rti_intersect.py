'''
Filter the results from bedtools intersect -wao between a set of exons and
retorcopy indels based on reciprocal overlap as a proportion of total
length. Output a file in the bedtools intersect output format.
'''

import sqlite3
import sys
import tempfile as tf
import subprocess

def get_intersect_dict(intersect_fname):
    column_names = [
        "exon_chromosome", "exon_start", "exon_end",
        "exon_name", "exon_score", "exon_strand",
        "rti_chromosome", "rti_start", "rti_end",
        "rti_name", "rti_score", "rti_strand",
        "overlap"
    ]
    intersect_dict = {}
    with open(intersect_fname, 'r') as f:
        for line in f:
            line_dict = {k:v for k,v in zip(column_names, line.strip().split())}

            if line_dict["rti_start"] == "-1":
                continue

            transcript_name = ".".join(line_dict["exon_name"].split(".")[0:-1])
            rti_name = line_dict["rti_name"]

            if transcript_name in intersect_dict:
                try:
                    intersect_dict[transcript_name][rti_name] \
                        += int(line_dict["overlap"])
                except KeyError:
                    intersect_dict[transcript_name][rti_name] \
                        = int(line_dict["overlap"])
            else:
                intersect_dict[transcript_name] \
                    = {rti_name:int(line_dict["overlap"])}

    return intersect_dict

def get_bed_dict(bed_fname):
    bed_dict = {}
    with open(bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            name = line_list[3]
            length = int(line_list[2]) - int(line_list[1])
            bed_dict[name] = {"length":length, "line":line.strip()}

    return bed_dict

def get_output_line(transcript, rti_name, exon_cur, rti_dict):
    transcript_list = list(
        exon_cur.execute(
            '''
            SELECT chromosome, start_coord, end_coord, transcript_id, strand
            FROM transcript WHERE transcript_id=?
            ''',
            (transcript, )).fetchone()
    )
    transcript_list.insert(-1,0)
    transcript_line = "\t".join(map(str, transcript_list))

    rti_line = rti_dict[rti_name]["line"]

    return (transcript_line, rti_line)

def filter_intersects(intersect_dict, exon_db, rti_bed, overlap_cutoff):
    '''
    For each transcript, check whether (total overlap/total exon length) is
    greater than the overlap cutoff. If it is, check the same for the
    RTI it overlaps with. If this also exceeds the overlap_cutoff,
    put the pair in the output list.
    '''

    out_intersects = []

    exon_conn = sqlite3.connect(exon_db)
    exon_cur = exon_conn.cursor()

    rti_dict = get_bed_dict(rti_bed)

    for transcript, intersects in intersect_dict.iteritems():
        try:
            total_exon_length = exon_cur.execute(
                "SELECT exon_length FROM transcript WHERE transcript_id=?",
                (transcript, )
            ).fetchone()[0]
        except:
            print(transcript)
            print(exon_cur.execute(
                "SELECT * FROM transcript WHERE transcript_id=?",
                (transcript, )
            ).fetchone())
            sys.exit()

        for rti_name, overlap in intersects.iteritems():
            rti_length = rti_dict[rti_name]["length"]

            exon_prop = float(overlap)/total_exon_length
            rti_prop = float(overlap)/rti_length

            if exon_prop >= overlap_cutoff and rti_prop >= overlap_cutoff:
                transcript_line, rti_line \
                    = get_output_line(transcript, rti_name, exon_cur, rti_dict)
                out_intersects.append([transcript_line, rti_line, str(overlap)])

    return out_intersects

def write_output_file(out_intersects, output_fname):
    with tf.NamedTemporaryFile() as tmp:
        for intersect in out_intersects:
            tmp.write("\t".join(intersect) + "\n")
        tmp.flush()

        subprocess.check_output(
            "sort -k1,1 -k2,2n {} > {}".format(tmp.name, output_fname),
            shell=True
        )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("intersect_fname")
    parser.add_argument("exon_db")
    parser.add_argument("rti_bed")
    parser.add_argument("overlap_cutoff", type=float)
    parser.add_argument("output_fname")
    args = parser.parse_args()

    intersect_dict = get_intersect_dict(args.intersect_fname)

    out_intersects = filter_intersects(
        intersect_dict, args.exon_db, args.rti_bed, args.overlap_cutoff)

    write_output_file(out_intersects, args.output_fname)
