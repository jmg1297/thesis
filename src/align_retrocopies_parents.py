'''
Inputs:
- mBED file of transcribed retroocpies with transcript, retrocopy, parent
- mBED file of untranscribed retrocopies with retrocopy, parent
- indexed FASTA file of retrocopy sequences
- indexed FASTA file of parent sequences
- output directory for alignment files
'''

import sys
import os
import sqlite3
import multiprocessing as mp
from Queue import Empty
import subprocess
import tempfile as tf
from math import floor
from random import shuffle

class Aligner(mp.Process):
    def __init__(self, queue, output_dir):
        mp.Process.__init__(self)
        self.queue = queue
        self.output_dir = output_dir
        self.align_cmd = "matcher {} {} {}"

    def run(self):
        while True:
            try:
                rc_name, retrocopy_fa_str, parent_fa_str = self.queue.get(timeout=10)
            except Empty:
                return 0
            else:
                self.align(rc_name, retrocopy_fa_str, parent_fa_str)

    def align(self, retrocopy_name, retrocopy_fa_str, parent_fa_str):
        output_path = os.path.join(self.output_dir, retrocopy_name + ".matcher")

        with tf.NamedTemporaryFile() as rc_fa, tf.NamedTemporaryFile() as parent_fa:
            rc_fa.write(retrocopy_fa_str)
            rc_fa.flush()
            parent_fa.write(parent_fa_str)
            parent_fa.flush()

            code = subprocess.check_output(
                "matcher {} {} {} 2>/dev/null".format(rc_fa.name, parent_fa.name, output_path),
                shell=True
            )


def get_idx_dict(idx_fname):
    idx_dict = {}
    with open(idx_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split(":")
            idx_dict[line_list[2]] = map(int, line_list[0:2])
    return idx_dict

def get_fa_str(lines, fasta):
    head = lines[1]
    tail = lines[1] - lines[0] + 1
    fa_str = subprocess.check_output(
        "head -{} {} | tail -{}".format(head, fasta, tail),
        shell=True
    )
    return fa_str

def make_alignments(
        retrocopy_parent_db, retrocopy_fasta, retrocopy_fasta_idx,
        parent_fasta, parent_fasta_idx, output_dir, num_procs, shuffle_pairs
    ):
    '''
    For each retrocopy, find its parent, get the corresponding sequences and put
    in temp files, then do a pairwise alignment and put the resulting file in
    the output directory.
    '''

    retrocopy_idx = get_idx_dict(retrocopy_fasta_idx)
    parent_idx = get_idx_dict(parent_fasta_idx)

    queue = mp.Queue()

    aligners = [Aligner(queue, output_dir) for _ in range(num_procs)]
    for a in aligners:
        a.start()

    #a = Aligner(queue, output_dir)

    conn = sqlite3.connect(retrocopy_parent_db)
    cur = conn.cursor()
    retrocopy_parents = cur.execute("SELECT * FROM retrogenes_parents").fetchall()

    if shuffle_pairs > 0:
        print("Shuffling pairs ... ")
        retrocopies = []
        parents = []
        for r in retrocopy_parents:
            if r is not None:
                retrocopies.append(r[0])
                parents.append(r[1])

        shuffle(parents)
        retrocopy_parents = [(r, p) for r, p in zip(retrocopies[0:shuffle_pairs], parents[0:shuffle_pairs])]


    conn.close()

    align_count = 0

    num_pairs = len(retrocopy_parents)
    percent = 0
    print('['+' '*50+']'),
    print('\r'),

    for i, (rc, parent) in enumerate(retrocopy_parents):

        if parent is None:
            continue

        retrocopy_fa_str = get_fa_str(retrocopy_idx[rc], retrocopy_fasta)

        try:
            pidx = parent_idx[parent]
        except KeyError:
            continue

        parent_fa_str = get_fa_str(pidx, parent_fasta)

        align_count += 1
        queue.put((rc, retrocopy_fa_str, parent_fa_str))

        #a.align(rc, retrocopy_fa_str, parent_fa_str)

        percent = 100.0*(float(i)/num_pairs)
        num_bars = int(floor(percent/2))
        print('\r'),
        print('['+'|'*num_bars+' '*(50-num_bars)+']'),
        sys.stdout.flush()

    print('\r'),
    print('['+'|'*50+']')
    sys.stdout.flush()

    print("Finished processing all pairs.")
    print("{} sent for alignment.".format(align_count))

    while(len(os.listdir(output_dir)) < align_count):
        print(
            "Waiting for {} alignment files"
                .format(align_count - len(os.listdir(output_dir)))
        )

    print("All alignments finished.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("retrocopy_parent_db")
    parser.add_argument("retrocopy_fasta")
    parser.add_argument("retrocopy_fasta_idx")
    parser.add_argument("parent_fasta")
    parser.add_argument("parent_fasta_idx")
    parser.add_argument("output_dir")
    parser.add_argument("num_procs", type=int)
    parser.add_argument("shuffle_pairs", type=int)
    args = parser.parse_args()

    make_alignments(
        args.retrocopy_parent_db,
        args.retrocopy_fasta, args.retrocopy_fasta_idx,
        args.parent_fasta, args.parent_fasta_idx,
        args.output_dir, args.num_procs, args.shuffle_pairs
    )
