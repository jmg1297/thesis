#!/usr/bin/python

import sqlite3
import os

def make_rsem_db(db_fname, samples_dir, filepath):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    cur.execute(
        '''
        CREATE TABLE stringtie_rsem (
        sample TEXT NOT NULL,
        transcript_id TEXT NOT NULL,
        gene_id TEXT NOT NULL,
        length INT NOT NULL,
        eff_length REAL NOT NULL,
        expected_count REAL NOT NULL,
        tpm REAL NOT NULL,
        fpkm REAL NOT NULL,
        isopct REAL NOT NULL
        )
        '''
    )

    for sample in os.listdir(samples_dir):
        fname = os.path.join(samples_dir, sample, filepath)
        data = []
        with open(fname, 'r') as f:
            num_cols = len(f.readline().strip().split())+1
            for line in f:
                row = line.strip().split()
                row[2] = int(row[2])
                row[3:] = map(float, row[3:])
                row = [sample] + row
                data.append(row)
        cur.executemany(
        '''
        INSERT INTO stringtie_rsem VALUES ({})
        '''.format(",".join(["?"]*num_cols)),
        data
        )

    conn.commit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("db_fname")
    parser.add_argument("samples_dir")
    parser.add_argument("filepath")
    args = parser.parse_args()

    make_rsem_db(args.db_fname, args.samples_dir, args.filepath)
