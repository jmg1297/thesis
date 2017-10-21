#!/usr/bin/python

'''
Make an SQLite DB of metadata for BLUEPRINT RNA-seq samples
'''

import sys
import json
import sqlite3

metadata_json = sys.argv[1]
metadata_db = sys.argv[2]

with open(metadata_json, 'r') as f:
    metadata = json.load(f)

conn = sqlite3.connect(metadata_db)
cur = conn.cursor()

try:
    cur.execute(
        '''
        CREATE TABLE metadata (
        code TEXT PRIMARY KEY,
        hybrid INT,
        strain TEXT,
        sex TEXT,
        cell TEXT
        )
        '''
    )
except sqlite3.OperationalError:
    cur.execute("DELETE FROM metadata")

hybrid_ints = {"pure":0, "hybrid":1}

data = []

for hybrid_status, strain_dict in metadata.iteritems():
    for strain, sex_dict in strain_dict.iteritems():
        for sex, cell_dict in sex_dict.iteritems():
            for cell, codes in cell_dict.iteritems():
                for code in codes:
                    data.append(
                        (code, hybrid_ints[hybrid_status], strain, sex, cell)
                    )

cur.executemany(
    "INSERT INTO metadata VALUES ({})".format(",".join(["?"]*len(data[0]))),
    data
)

conn.commit()
