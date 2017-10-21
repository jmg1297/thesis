'''
Create a DB from a BED file of retrocopies and a BED file of their constituent
blocks, with one table for each file.
'''

import sqlite3
from collections import OrderedDict

def make_retrocopy_block_db(retrocopy_fname, block_fname, db_fname):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    cur.execute(
        '''
        CREATE TABLE retrocopy (
            chromosome text,
            start_coord integer,
            end_coord integer,
            name text,
            strand text,
            num_blocks integer,
            block_length integer,
            PRIMARY KEY (name)
        )
        '''
    )

    cur.execute(
        '''
        CREATE TABLE block (
            chromosome text,
            start_coord integer,
            end_coord integer,
            name text,
            block_number integer,
            strand text,
            PRIMARY KEY (name, block_number)
        )
        '''
    )

    retrocopy_dict = OrderedDict()
    block_data = []

    with open(retrocopy_fname, 'r') as retrocopies:
        for line in retrocopies:
            line_list = line.strip().split()
            name = line_list[3]
            retrocopy_dict[name] = OrderedDict([
                ("chromosome", line_list[0]),
                ("start_coord", int(line_list[1])),
                ("end_coord", int(line_list[2])),
                ("name", name),
                ("strand", line_list[5]),
                ("num_blocks", 0),
                ("block_length", 0)
            ])

    with open(block_fname, 'r') as blocks:
        for line in blocks:
            line_list = line.strip().split()
            chromosome = line_list[0]
            start_coord = int(line_list[1])
            end_coord = int(line_list[2])
            name, block_id = line_list[3].split(":")[2:]
            block_number = int(block_id.split(".")[1])
            strand = line_list[5]
            block_data.append(
                (chromosome, start_coord, end_coord, name, block_number, strand)
            )
            retrocopy_dict[name]["num_blocks"] += 1
            retrocopy_dict[name]["block_length"] += (end_coord - start_coord)

    retrocopy_data = [tuple(d.values()) for d in retrocopy_dict.values()]

    retrocopy_cols = len(retrocopy_data[0])
    block_cols = len(block_data[0])

    cur.executemany(
        "INSERT INTO retrocopy VALUES ({})"\
            .format(",".join(["?"]*retrocopy_cols)),
        retrocopy_data
    )

    cur.executemany(
        "INSERT INTO block VALUES ({})"\
            .format(",".join(["?"]*block_cols)),
        block_data
    )

    conn.commit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("retrocopy_fname")
    parser.add_argument("block_fname")
    parser.add_argument("db_fname")
    args = parser.parse_args()

    make_retrocopy_block_db(
        args.retrocopy_fname, args.block_fname, args.db_fname)
