'''
Given the output of bedtools merge applied to matches of mm10 retrocopies in
CAST, split the merged regions based on whether the retrocopy parents agree
'''

import logging
import sqlite3

logging.basicConfig(level=logging.INFO)

def split_by_parent_agreement(merge_fname, rc_parent_db, agree_out_fname, disagree_out_fname, none_out_fname):

    conn = sqlite3.connect(rc_parent_db)
    cur = conn.cursor()

    get_parent_query = '''
        SELECT parents FROM retrogenes_parents WHERE retrogene IN ({})
    '''

    logging.info("Starting to check merged matches")

    with open(merge_fname, 'r') as f,\
            open(agree_out_fname, 'wa') as aout,\
            open(disagree_out_fname, 'wa') as dout,\
            open(none_out_fname, 'wa') as nout:
        for i,line in enumerate(f):
            if (i+1)%1000 == 0:
                logging.info("{} lines processed".format(i+1))

            retrocopies = line.strip().split()[4].split(",")
            rows = cur.execute(get_parent_query.format(",".join(["?"]*len(retrocopies))), retrocopies).fetchall()
            parents = sorted(list(set([r[0] for r in rows])))
            out_line = line.strip() + "\t" + ",".join(map(str, parents)) + "\n"
            if None in parents:
                nout.write(out_line)
            elif len(parents) == 1:
                aout.write(out_line)
            else:
                dout.write(out_line)
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("merge_fname")
    parser.add_argument("rc_parent_db")
    parser.add_argument("agree_out_fname")
    parser.add_argument("disagree_out_fname")
    parser.add_argument("none_out_fname")
    args = parser.parse_args()

    split_by_parent_agreement(
        args.merge_fname,
        args.rc_parent_db,
        args.agree_out_fname,
        args.disagree_out_fname,
        args.none_out_fname
    )
