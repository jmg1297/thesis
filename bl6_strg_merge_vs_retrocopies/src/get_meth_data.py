'''
Given a list of BED files with labels and a list of BEDgraph files of
methylation data, do the following:
for each BED file:
  for each region:
    for each BEDgraph file:
      extract region's methylation levels
      calculate summary statistics (overall level, mean, std dev)
      return None if no CpGs in region
Store results in a dict and write to JSON
'''

import sys
import json
import os
import linecache
from math import *
import subprocess
import multiprocessing as mp
from Queue import Empty
import logging

import numpy as np

SCRIPTDIR="/home/jg600/bs-seq.bedgraphs/scripts"
logging.basicConfig(level=logging.INFO)

class GetMethProcess(mp.Process):
    def __init__(self, job_queue, results_queue):
        mp.Process.__init__(self)
        self.job_queue = job_queue
        self.results_queue = results_queue
        logging.info("Initialised {}".format(self.name))

    def run(self):
        logging.info("{} starting to retrieve jobs from queue".format(self.name))

        while True:
            job = self.job_queue.get()
            if job == "STOP":
                logging.info("Received STOP signal, stopping {}".format(self.name))
                return 0
            else:
                l, bed, gl, bedgraph_fname, idx_dict = job
                logging.info("Retrieved {}/{} pair, processing".format(l, gl))
                data = self.get_meth_data(bed, bedgraph_fname, idx_dict)
                logging.info("Finished {}/{} pair, putting results on queue".format(l, gl))
                self.results_queue.put([l, gl, data])


    def get_meth_data(self, bed, bedgraph_fname, idx_dict):
        data = {}
        num_regions = len(bed)
        for i, (key, info) in enumerate(bed.iteritems()):
            if (i+1)%100 == 0:
                logging.info("{}: Processed {}/{} regions".format(self.name, i+1, num_regions))

            data[key] = {}
            chrom = info["chrom"]
            region_start = info["start"]
            region_end = info["end"]

            first_line, last_line = self.bisection_search(
                bedgraph_fname,
                idx_dict[chrom]["first"],
                idx_dict[chrom]["last"],
                region_start,
                region_end
            )

            if first_line is None:
                data[key] = None
            else:
                meth_stats = self.get_meth_stats(
                    bedgraph_fname, first_line, last_line)
                if None in meth_stats:
                    data[key] = None
                else:
                    data[key] = {
                        "mean":meth_stats[0],
                        "variance":meth_stats[1],
                        "level":meth_stats[2]
                    }

        logging.info("{}: Processed {}/{} regions".format(self.name, num_regions, num_regions))
        return data

    def bisection_search(self, bedgraph_fname, chrom_first_line, chrom_last_line, region_start, region_end):
        '''
        Search through a bedgraph to find the first and last lines contained within
        a given genomic region.
        '''

        # Get the data from the first and last lines of the region to be searched
        first_line_list = linecache.getline(bedgraph_fname, chrom_first_line).strip().split()
        last_line_list = linecache.getline(bedgraph_fname, chrom_last_line).strip().split()

        region_first_line = 0
        region_last_line = 0

        if region_start <= int(first_line_list[1]):
            # Region starts before first CpG in the chromosome
            region_first_line = chrom_first_line
        else:
            found = False
            top_line = chrom_first_line
            bottom_line = chrom_last_line
            mid_line = int(floor((top_line + bottom_line)/2))

            while (bottom_line != top_line + 1) and not found:
                line = linecache.getline(bedgraph_fname, mid_line).strip().split()
                position = int(line[1])

                if position > region_start:
                    # Take top half
                    bottom_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position < region_start:
                    # Take bottom half
                    top_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position == region_start:
                    # Exact hit
                    region_first_line = mid_line
                    found = True

            if not found:
                region_first_line = bottom_line

        if region_end >= int(last_line_list[1]):
            # Region ends after last CpG in the chromosome
            region_last_line = chrom_last_line
        else:
            found = False
            top_line = region_first_line
            bottom_line = chrom_last_line
            mid_line = int(floor((top_line + bottom_line)/2))

            while (bottom_line != top_line + 1) and not found:
                line = linecache.getline(bedgraph_fname, mid_line).strip().split()
                position = int(line[1])

                if position > region_end:
                    # Take top half
                    bottom_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position < region_end:
                    # Take bottom half
                    top_line = mid_line
                    mid_line = int(floor((top_line + bottom_line)/2))
                elif position == region_end:
                    region_last_line = mid_line
                    found = True

            if not found:
                region_last_line = top_line

        if region_first_line == region_last_line:
            # Either there are no CpGs in this region, or this line contains the
            # only CpG in the region
            line = linecache.getline(bedgraph_fname, region_first_line).strip().split()
            position = int(line[1])
            if region_start <= position <= region_end:
                # Only CpG in the region
                pass
            else:
                # No CpGs in the region
                return None,None

        return [region_first_line, region_last_line]

    def get_meth_stats(self, bedgraph_fname, first_line, last_line):
        '''
        Given a pair of lines in a bedgraph file, collect all the data in between
        these lines and calculate summary statistics
        '''
        meth_levels = []
        meth_reads = []
        cvg = []

        for i in range(first_line, last_line + 1):
            line_list = linecache.getline(bedgraph_fname, i).strip().split()
            try:
                meth_levels.append(float(line_list[3]))
                mr = int(line_list[4])
                umr = int(line_list[5])
                meth_reads.append(mr)
                cvg.append(mr + umr)
            except ValueError:
                # No data for this CpG
                continue

        if meth_levels == []:
            return None, None, None

        mean = np.mean(meth_levels)
        variance = np.std(meth_levels)**2
        level = float(sum(meth_reads))/sum(cvg)

        return mean, variance, level

def open_bed(bed_fname):
    logging.info("Reading in {}".format(bed_fname))

    bed = {}
    with open(bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            key = ":".join(line_list[0:4])
            bed[key] = {
                "chrom":line_list[0],
                "start":int(line_list[1]),
                "end":int(line_list[2])
            }

    logging.info("Finished reading {}".format(bed_fname))
    return bed

def get_bedgraph_idx(idx_fname):
    logging.info("Fetching bedgraph index info from {}".format(idx_fname))

    index_dict = {}
    with open(idx_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split(",")
            index_dict[line_list[0]] = {
                "first":int(line_list[1]),
                "last":int(line_list[2])
            }

    logging.info("Finished")
    return index_dict

def get_meth_data(input_beds, bed_labels, bedgraphs, bedgraph_labels, num_procs):
    logging.info("Starting to get methylation data")

    meth_data = {l:{gl:{} for gl in bedgraph_labels} for l in bed_labels}
    job_queue = mp.Queue()
    results_queue = mp.Queue()

    logging.info("Starting data fetching processes")
    procs = [GetMethProcess(job_queue, results_queue) for _ in range(num_procs)]
    for p in procs:
        p.start()

    logging.info("All processes started")

    logging.info("Checking bedgraph index files")

    bedgraph_idx_fnames = []
    for graph_fname in bedgraphs:
        graph_idx_fn = graph_fname + ".idx"
        if not os.path.isfile(graph_idx_fn):
            logging.info("Missing index file for {}".format(graph_fname))
            logging.info("Creating index file now ... ")
            subprocess.call(
                "bash {}/index_bedgraph.sh {}".format(SCRIPTDIR, graph_fname),
                shell=True
            )
            logging.info("Finished making index file")
        bedgraph_idx_fnames.append(graph_idx_fn)

    logging.info("Checked all bedgraph index files")

    logging.info("Starting to put jobs on queue")
    num_results = len(input_beds)*len(bedgraphs)

    for bed_fname, l in zip(input_beds, bed_labels):
        bed = open_bed(bed_fname)
        for graph_fname, gl, idx_fname in zip(bedgraphs, bedgraph_labels, bedgraph_idx_fnames):
            idx_dict = get_bedgraph_idx(idx_fname)

            logging.info("Putting {}/{} pair on queue".format(l, gl))
            job_queue.put([l, bed, gl, graph_fname, idx_dict])

    for _ in range(num_procs):
        job_queue.put("STOP")

    logging.info("All jobs on queue")
    logging.info("Retrieving results")

    count = 0
    while True:
        l, gl, data = results_queue.get()
        count += 1
        meth_data[l][gl] = data
        if count == num_results:
            break

    logging.info("All results retrieved")
    return meth_data

def write_data(meth_data, output_fname):
    logging.info("Writing methylation data to {}".format(output_fname))
    with open(output_fname, 'wa') as out:
        json.dump(meth_data, out, indent=4, sort_keys=True)
    logging.info("Finished writing data")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputBed", action="append")
    parser.add_argument("-l", "--bedLabel", action="append")
    parser.add_argument("-g", "--BEDgraph", action="append")
    parser.add_argument("-L", "--graphLabel", action="append")
    parser.add_argument("num_procs", type=int)
    parser.add_argument("output_fname")
    args = parser.parse_args()

    meth_data = get_meth_data(
        args.inputBed, args.bedLabel,
        args.BEDgraph, args.graphLabel,
        args.num_procs
    )

    write_data(meth_data, args.output_fname)
