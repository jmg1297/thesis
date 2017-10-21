'''
Create plots showing histone modification peaks at various distances from
expressed RTIs compared to randomly chosen RTIs. For each histone modification,
create a subplot. Each subplot is as follows. For each distance from the RTI,
draw a boxplot of the distribution from the random data, and a cross for the
data from the expressed RTIs.
'''

import matplotlib.pyplot as plt
import multiprocessing as mp
from Queue import Empty
import subprocess
import tempfile as tf
import logging

logging.basicConfig(level=logging.INFO)

class IntersectProcess(mp.Process):
    def __init__(self, job_queue, results_queue, close_queue=None):
        mp.Process.__init__(self)
        self.job_queue = job_queue
        self.results_queue = results_queue
        self.close_queue = close_queue
        self.intersect_cmd = "bedtools window -u -w {} -a {} -b {} | wc -l"

    def run(self):
        while True:
            try:
                distance, rti_fn, hmod_fn, hmod_label = self.job_queue.get(timeout=10)
            except Empty:
                return 0
            else:
                self.intersect(distance, rti_fn, hmod_fn, hmod_label)

    def intersect(self, distance, rti_fn, hmod_fn, hmod_label):

        logging.debug(self.intersect_cmd.format(distance, rti_fn, hmod_fn))
        num_intersects = int(
            subprocess.check_output(
                self.intersect_cmd.format(distance, rti_fn, hmod_fn),
                shell=True
            )
        )

        if self.close_queue is not None:
            self.close_queue.put(rti_fn)

        self.results_queue.put((hmod_label, distance, num_intersects))

def get_random_data(
        all_rti_bed, sample_size, num_samples, hmod_bed_list,
        hmod_label_list, distance_list, num_procs):
    '''
    Get intersection data for randomly chosen subsets of all of the RTIs
    '''
    logging.info("Getting random data")

    hmod_dict = {h:fn for h,fn in zip(hmod_label_list, hmod_bed_list)}

    random_data = {h:{d:[] for d in distance_list} for h in hmod_label_list}

    job_queue = mp.Queue()
    results_queue = mp.Queue()
    # Used so IntersectProcess can tell main to close tmpfiles, otherwise
    # program crashes because too many files are open
    close_queue = mp.Queue()

    logging.info("Spawning intersection processes")

    procs = [
        IntersectProcess(job_queue, results_queue, close_queue=close_queue)
        for _ in range(num_procs)
    ]
    
    for p in procs:
        p.start()

    logging.info("All intersection processes running")

    tmpfiles = []
    tmpdict = {}

    for i in range(num_samples):
        if (i)%100 == 0:
            logging.info("Processed {} random samples".format(i))
        for h,hfn in hmod_dict.iteritems():
            for d in distance_list:
                tmpfiles.append(tf.NamedTemporaryFile())
                tmpdict[tmpfiles[-1].name] = tmpfiles[-1]

                subprocess.call(
                    "shuf -n {} {} | sort -k1,1 -k2,2n > {}"
                        .format(sample_size, all_rti_bed, tmpfiles[-1].name),
                    shell=True
                )
                tmpfiles[-1].flush()
                job_queue.put((d, tmpfiles[-1].name, hfn, h))
                try:
                    tmpname = close_queue.get(timeout=0.05)
                except Empty:
                    pass
                else:
                    tmpdict[tmpname].close()

    logging.info("All random samples on job queue")

    expected_results = len(distance_list)*len(hmod_label_list)*num_samples

    logging.info("Waiting for {} results".format(expected_results))

    results_count = 0
    while results_count < expected_results:
        h,d,n = results_queue.get()
        random_data[h][d].append(n)
        results_count += 1
        if results_count%100 == 0:
            logging.info(
                "Retrieved {} results from random samples".format(results_count)
            )

    logging.info("Retrieved all results")

    logging.info("Returning random data")
    return random_data

def get_expr_data(
        expr_rti_bed, hmod_bed_list, hmod_label_list, distance_list, num_procs):
    '''
    Get intersection data for expressed RTIs
    '''
    logging.info("Getting data from expressed RTIs")
    hmod_dict = {h:fn for h,fn in zip(hmod_label_list, hmod_bed_list)}

    expr_data = {h:{d:0 for d in distance_list} for h in hmod_label_list}

    job_queue = mp.Queue()
    results_queue = mp.Queue()

    logging.info("Spawning intersection processes")

    procs = [IntersectProcess(job_queue, results_queue) for _ in range(num_procs)]
    for p in procs:
        p.start()

    logging.info("All intersection processes running")

    for h,hfn in hmod_dict.iteritems():
        for d in distance_list:
            job_queue.put((d, expr_rti_bed, hfn, h))

    logging.info("All combinations on job queue")

    expected_results = len(hmod_label_list)*len(distance_list)

    logging.info("Waiting for {} results".format(expected_results))

    results_count = 0
    while results_count < expected_results:
        h,d,n = results_queue.get()
        expr_data[h][d] = n
        results_count += 1

    logging.info("Retrieved all results")

    logging.info("Returning expressed RTI data")
    return expr_data

def plot_intersects(
        expr_data, random_data, hmod_label_list, distance_list,
        plot_title, img_fname):
    '''
    For each histone modification, draw a subplot. Within each subplot, for each
    distance draw a boxplot of the random data and add a cross for the
    expressed data.
    '''

    logging.info("Drawing plots")

    fig, axarr = plt.subplots(1, len(hmod_label_list), sharey=True)
    fig.set_size_inches(5*len(hmod_label_list), 6)
    if len(hmod_label_list) == 1:
        axarr = [axarr]

    ax_dict = {h:ax for h,ax in zip(hmod_label_list, axarr)}

    plt.suptitle(plot_title)

    distance_list.sort()

    for h in hmod_label_list:
        logging.info("Drawing plot for {}".format(h))
        ax = ax_dict[h]
        subplot_data = random_data[h]
        data = [subplot_data[d] for d in distance_list]
        ax.boxplot(data, labels=distance_list)
        ax.set_xticklabels(distance_list, rotation=60)
        ax.set_title(h)

        expr_points = expr_data[h]
        for i,d in enumerate(distance_list):
            ax.plot(i+1, expr_points[d], 'rx', ms=10, linewidth=0)

    axarr[0].set_ylabel("Number of RTIs")
    axarr[0].set_ylim([0,300])

    logging.info("Finished all plots, saving figure")

    plt.savefig(img_fname, bbox_inches="tight")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("expr_rti_bed")
    parser.add_argument("all_rti_bed")
    parser.add_argument("sample_size", type=int)
    parser.add_argument("num_samples", type=int)
    parser.add_argument("num_procs", type=int)
    parser.add_argument("-b", "--hmod_bed", action="append")
    parser.add_argument("-l", "--hmod_label", action="append")
    parser.add_argument("-d", "--distance", action="append", type=int)
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    random_data = get_random_data(
        args.all_rti_bed,
        args.sample_size,
        args.num_samples,
        args.hmod_bed,
        args.hmod_label,
        args.distance,
        args.num_procs
    )

    expr_data = get_expr_data(
        args.expr_rti_bed,
        args.hmod_bed,
        args.hmod_label,
        args.distance,
        args.num_procs
    )

    plot_intersects(
        expr_data,
        random_data,
        args.hmod_label,
        args.distance,
        args.plot_title,
        args.img_fname
    )
