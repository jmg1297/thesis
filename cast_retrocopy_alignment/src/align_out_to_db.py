'''
Create an SQLite database from the results of a blastn alignment
'''

import sqlite3

def make_row(line_list):
    row = []
    for i,x in enumerate(line_list):
        if i in [0,1]:
            row.append(x)
        elif i in [2,10,11]:
            row.append(float(x))
        else:
            row.append(int(x))
    return row

def make_db(align_out_fname, db_fname):

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    cur.execute(
        '''
        CREATE TABLE alignment (
            query_name TEXT,
            subject_name TEXT,
            percent_identity REAL,
            alignment_length INT,
            mismatches INT,
            gap_opens INT,
            query_start INT,
            query_end INT,
            subject_start INT,
            subject_end INT,
            expect_value REAL,
            bit_score REAL
        )
        '''
    )

    data = []
    with open(align_out_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            data.append(make_row(line_list))

    cur.executemany(
        "INSERT INTO alignment VALUES ({})".format(",".join(["?"]*len(data[0]))),
        data
    )

    conn.commit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("align_out_fname")
    parser.add_argument("db_fname")
    args = parser.parse_args()

    make_db(args.align_out_fname, args.db_fname)
