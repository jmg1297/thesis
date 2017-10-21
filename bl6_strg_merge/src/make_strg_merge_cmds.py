#!/usr/bin/python

import os
import sqlite3

def get_gtf_lists(metadata_db, strain, hybrid, strg_dir, gtf_fname):
    conn = sqlite3.connect(metadata_db)
    cur = conn.cursor()

    sample_lists = {}

    query = '''
    SELECT DISTINCT(cell) FROM metadata WHERE hybrid=? AND strain=?
    '''
    cell_types = [r[0] for r in cur.execute(query, (hybrid, strain)).fetchall()]

    query = '''
    SELECT DISTINCT(sex) FROM metadata WHERE hybrid=? AND strain=?
    '''
    sexes = [r[0] for r in cur.execute(query, (hybrid, strain)).fetchall()]

    for cell in cell_types:
        query = '''
        SELECT code, strain, sex, cell FROM metadata
        WHERE hybrid=? AND strain=? AND cell=?
        '''
        sample_lists[cell] = ["_".join(r) for r in cur.execute(query, (hybrid, strain, cell)).fetchall()]

    for sex in sexes:
        query = '''
        SELECT code, strain, sex, cell FROM metadata
        WHERE hybrid=? AND strain=? AND sex=?
        '''
        sample_lists[sex] = ["_".join(r) for r in cur.execute(query, (hybrid, strain, sex)).fetchall()]

    for sex in sexes:
        for cell in cell_types:
            query = '''
            SELECT code, strain, sex, cell FROM metadata
            WHERE hybrid=? AND strain=? AND sex=? AND cell=?
            '''
            key = sex + "_" + cell
            sample_lists[key] = ["_".join(r) for r in cur.execute(query, (hybrid, strain, sex, cell)).fetchall()]

    query = '''
    SELECT code, strain, sex, cell FROM metadata
    WHERE hybrid=? AND strain=?
    '''
    sample_lists["ALL"] = ["_".join(r) for r in cur.execute(query, (hybrid, strain)).fetchall()]

    gtf_lists = {
        key:[os.path.join(strg_dir, sample, gtf_fname) for sample in samples]
        for key, samples in sample_lists.iteritems()
    }

    return gtf_lists

def write_cmds(gtf_lists, merge_outdir, output_file):

    cmd_str = "stringtie --merge -o {} -m 200 -l {}_MSTRG {}"

    with open(output_file, 'wa') as f:
        for key, gtfs in gtf_lists.iteritems():
            out_path = os.path.join(merge_outdir, "{}_merge.gtf".format(key))
            f.write(cmd_str.format(out_path, key, " ".join(gtfs)) + "\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_db")
    parser.add_argument("strain")
    parser.add_argument("hybrid", type=int)
    parser.add_argument("strg_dir")
    parser.add_argument("gtf_fname")
    parser.add_argument("merge_outdir")
    parser.add_argument("output_file")
    args = parser.parse_args()

    gtf_lists = get_gtf_lists(
        args.metadata_db, args.strain, args.hybrid,
        args.strg_dir, args.gtf_fname
    )

    write_cmds(gtf_lists, args.merge_outdir, args.output_file)
