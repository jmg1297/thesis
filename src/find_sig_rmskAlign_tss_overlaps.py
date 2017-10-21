'''
Script to identify rmskAlignBaseline elements that are significantly enriched
for RTI transcript TSSs by generating a random null distribution and comparing
with the actual values. The null distribution is generated as follows:
    - Input: intersection of consensus RTI transcripts and RTIs
    - From this, get the distributions of RTI content classification, transcript
      overhang/underhang, and transcript strand
    - Input: all RTIs with content classification
    - Pick random subsets of RTIs, following the distribution of content
      classification.
    - Generate random TSSs for these RTIs using the distributions of overhang/
      underhang and strand
    - Intersect the random TSSs with the rmskAlign coordinates and extract
      the number of appearances for each type of rmskAlign element
    - Repeat as many times as necessary
    - Input: intersection of actual RTI TSSs and rmskAlign
    - Get the numbers of appearances for each type of rmskAlign element from the
      intersection with the actual TSSs, filtering out those representing less
      than x% of the total.
    - For the types of rmskAlign element remaining, fetch the null distributions,
      calculate the mean and SD for each, and normalise to a (0,1) distribution.
    - Normalise the actual value and calculate the number of SDs it lies away
      from the mean
    - Report those lying more than p standard deviations from 0.
    - Visualise results in a boxplot or similar.
'''

import sys
import os
import subprocess
from StringIO import StringIO
import re
import time
import multiprocessing as mp
import random
import numpy.random as nprand
import scipy.stats as spstats
import tempfile as tf
import matplotlib.pyplot as plt
import numpy as np
import json
from math import ceil

def get_flike(cmd_str):
    flike = StringIO(subprocess.check_output(cmd_str, shell=True))
    return flike

def get_props(flike, total):
    counts_dict = {}
    for line in flike:
        key, count = line.strip().split()
        counts_dict[key] = float(count)

    props_dict = {k:float(count)/total \
                    for k,count in counts_dict.iteritems()}
    return props_dict

def get_transcript_feat_distns(consensus_rti_intersect_fname):

    print("Getting transcript feature distributions ... ")
    sys.stdout.flush()

    print("Getting total number of RTIs ... "),
    sys.stdout.flush()

    total_rtis = int(
        subprocess.check_output(
            '''
            cat %s | awk '{print $7,$8,$9,$10}' | sort | uniq | wc -l
            ''' % consensus_rti_intersect_fname,
            shell=True
        )
    )

    print("Done")
    sys.stdout.flush()

    print("Getting total number of transcripts ... "),
    sys.stdout.flush()

    total_transcripts = int(
        subprocess.check_output(
            '''
            cat %s | awk '{print $1,$2,$3,$4}' | sort | uniq | wc -l
            ''' % consensus_rti_intersect_fname,
            shell=True
        )
    )

    print("Done")
    sys.stdout.flush()

    print("Getting RTI classification distribution ... "),
    sys.stdout.flush()

    flike = get_flike(
            '''
            cat %s | awk '{print $13}' | sort | uniq -c | awk '{print $2,$1}'
            ''' % consensus_rti_intersect_fname
    )

    classification_props = get_props(flike, total_rtis)

    print("Done")
    print("Getting strand distribution ... "),
    sys.stdout.flush()

    flike = get_flike(
            '''
            cat %s | awk '{print $1,$2,$3,$4,$5,$6}' | sort | uniq |
            awk '{print $6}' | sort | uniq -c | awk '{print $2,$1}'
            ''' % consensus_rti_intersect_fname
    )

    strand_props = get_props(flike, total_transcripts)

    print("Done")
    print("Getting transcript relative positions ... "),
    sys.stdout.flush()

    flike = get_flike(
            "cat %s | awk '{print $3-$2, $2-$8}'" % consensus_rti_intersect_fname
    )

    relative_positions = [tuple(map(int, line.strip().split())) \
                            for line in flike]

    print("Done")
    print("All transcript feature distributions extracted")
    sys.stdout.flush()

    return classification_props, strand_props, relative_positions

def get_rtis_by_classification(content_classified_rtis_fname):

    print("Getting dict of RTIs by classification ... "),
    sys.stdout.flush()

    rtis_by_classification = {}

    with open(content_classified_rtis_fname, 'r') as f:
        for line in f:
            line_list = re.split("\s+", line.strip())
            classification = line_list[6]
            classification = re.sub("\?", "", classification)
            try:
                rtis_by_classification[classification].append(line.strip())
            except KeyError:
                rtis_by_classification[classification] = [line.strip()]

    print("Done")
    sys.stdout.flush()

    return rtis_by_classification

class NullGenProcess(mp.Process):

    input_queue = None
    summary_queue = None
    rtis_by_classification = None
    classification_sample_counts = None
    strand_props = None
    relative_positions = None
    rmsk_align_bed_fname = None
    min_prop = None
    element_class_dict = None

    def __init__(self):
        mp.Process.__init__(self)
        self.daemon = True

    def run(self):
        while True:
            token = self.input_queue.get()
            if token == "STOP":
                return None
            else:
                sample_summary = self.generate_sample()
                self.summary_queue.put(sample_summary)

    def generate_sample(self):
        rti_sample = self.get_rti_sample()
        random_tss_sample = self.generate_rand_tss(rti_sample)
        intersect_summary = self.intersect(random_tss_sample)
        return intersect_summary

    def get_rti_sample(self):
        rti_sample = []
        for classification, count in self.classification_sample_counts.iteritems():
            rti_sample += random.sample(
                self.rtis_by_classification[classification], count
            )

        return rti_sample

    def jitter(self, n):
        return n + random.randint(-10, 10)

    def generate_rand_tss(self, rti_sample):
        random_tss_sample = []

        strand_opts = []
        strand_p = []

        for k,v in self.strand_props.iteritems():
            strand_opts.append(k)
            strand_p.append(v)

        for line in rti_sample:
            line_list = re.split("\s+", line)
            chrom = line_list[0]
            rel_start, length = random.choice(self.relative_positions)
            rel_start = self.jitter(rel_start)
            length = self.jitter(length)
            start = int(line_list[1]) + rel_start
            end = start + length
            strand = nprand.choice(strand_opts, p=strand_p)
            if strand == ".":
                strand = random.choice(["+", "-"])

            if strand == "+":
                tss = str(start)
            elif strand == "-":
                tss = str(end)

            random_tss_sample.append("\t".join([chrom, tss, tss]))

        return random_tss_sample

    def intersect(self, random_tss_sample):
        tmp = tf.NamedTemporaryFile()
        tmp.write("\n".join(random_tss_sample) + "\n")
        tmp.flush()

        intersect_cmd = "bedtools intersect -wao -a %s -b %s" \
                            % (tmp.name, self.rmsk_align_bed_fname)

        intersect_output = subprocess.check_output(intersect_cmd, shell=True)

        tmp.close()

        tmp_intersect = tf.NamedTemporaryFile()
        tmp_intersect.write(intersect_output)
        tmp_intersect.flush()

        intersect_summary = get_intersect_summary(tmp_intersect.name, 7, self.min_prop, self.element_class_dict)

        tmp_intersect.close()

        return intersect_summary

def combine_summaries(summaries):
    null_distns = {}
    for summ in summaries:
        for element, count in summ.iteritems():
            try:
                null_distns[element].append(count)
            except KeyError:
                null_distns[element] = [count]
    return null_distns

def generate_null_distns(
        transcript_feat_distns, rtis_by_classification, rmsk_align_bed_fname,
        sample_size, num_samples, min_prop, num_procs, element_class_dict):

    print("Generating null distributions ... ")
    sys.stdout.flush()

    classification_props, strand_props, relative_positions = transcript_feat_distns

    classification_sample_counts = {k:int(classification_props[k]*sample_size) \
                                        for k in classification_props.keys()}

    input_queue = mp.Queue()
    summary_queue = mp.Queue()

    for name, value in locals().iteritems():
        setattr(NullGenProcess, name, value)

    gen_procs = [NullGenProcess() for _ in range(num_procs)]
    for p in gen_procs:
        p.start()

    summaries = []

    for i in (range(num_samples) + ["STOP"]*num_procs):
        input_queue.put(i)

    while len(summaries) < num_samples:
        print("{}/{} summaries received\r".format(len(summaries), num_samples)),
        sys.stdout.flush()
        summaries.append(summary_queue.get())

    print("{}/{} summaries received\r".format(len(summaries), num_samples))

    for p in gen_procs:
        p.terminate()

    null_distns = combine_summaries(summaries)

    print("Finished null distributions")
    return null_distns

def get_intersect_summary(intersect_fname, element_field, min_prop, element_class_dict):
    summary = {}
    flike = StringIO(
        subprocess.check_output(
            '''
            cat %s | awk '{print $%d}' | cut -d: -f1,2,3 | sort | uniq -c
            ''' % (intersect_fname, int(element_field)),
            shell=True
        )
    )

    for line in flike:
        count, element = line.strip().split()
        try:
            element_class = element_class_dict[element]
        except KeyError:
            element_class = element

        try:
            summary[element_class] += int(count)
        except KeyError:
            summary[element_class] = int(count)

    filtered_summary = {}

    total = sum(summary.values())
    for element,count in summary.iteritems():
        if float(count)/total >= min_prop:
            filtered_summary[element] = count

    return filtered_summary

def get_params(vals):
    return {"mean":float(np.mean(vals)), "sd":float(np.std(vals))}

def standardise(vals, params):
    #vals = np.array(vals)
    #standard = (vals - params["mean"])/params["sd"]
    #return standard.tolist()
    return vals

def get_p_value(query_data, standard_null_distn):
    sorted_null = sorted(standard_null_distn, reverse=True)
    for i,v in enumerate(sorted_null):
        if query_data > v:
            return float(i)/len(sorted_null)
    return 1.0

def draw_plot(
        standard_null_distns, standard_query_values, query_distance,
        null_distn_params, sample_size, sig_dist, img_fname):

    # boxplots, one for each element
    # order boxplots L to R according to distance of query from mean
    # mark query value with cross
    # horizontal line to show sig_dist
    # p value above significant ones

    elements = sorted(query_distance.keys(), key=lambda k: null_distn_params[k]["mean"])

    fig, ax = plt.subplots(figsize=(ceil(len(elements)/4.0), 8))
    #fig, ax = plt.subplots(figsize=(10,4))
    #fig, ax = plt.subplots()
    plt.title("RTI Transcript TSSs vs Repeat Elements")

    if "." in elements:
        elements.remove(".")
    
    #hist_null_distns(standard_null_distns, elements)

    random_data = [standard_null_distns[e] for e in elements] # float
    query_data = [standard_query_values[e] for e in elements] # float

    ax.boxplot(random_data, labels = elements)

    for i,v in enumerate(query_data):
        ax.plot(i+1, v, 'rx', ms=10, linewidth=0)

    for i,e in enumerate(elements):
        p = get_p_value(query_data[i], standard_null_distns[e])
        if p <= 0.05:
        #if query_distance[e] >= sig_dist:
            #p = float(spstats.norm.sf(query_distance[e]))
            ax.text(
                i+1, query_data[i]+np.sign(query_data[i])*1.5, "$p=%.3f$" % p,
                #i+1, query_distance[e] + 1.5, "$p=%.3f$" % p,
                horizontalalignment='center', fontsize=6
            )

    xmin, xmax = ax.get_xlim()
    #ax.plot((xmin, xmax), (sig_dist, sig_dist), 'k--')
    #ax.plot((xmin, xmax), (-sig_dist, -sig_dist), 'k--')

    ymin, ymax = ax.get_ylim()
    ax.set_ylim((ymin-1, ymax+1))

    ax.set_xticklabels(elements, rotation=90)

    ax.set_xlabel("RepeatMasker element")
    ax.set_ylabel("Normalised number of TSSs")

    plt.savefig(img_fname, bbox_inches='tight', pad_inches=5)

def hist_null_distns(standard_null_distns, elements):

    try:
        os.mkdir("null_hists")
    except:
        pass

    for element in elements:
        data = standard_null_distns[element]
        fig, ax = plt.subplots()
        ax.hist(data, lw=0)

        plt.savefig("null_hists/%s.svg" % element)

        plt.close(fig)

def plot_significant_elements(
        query_values, null_distns,
        sample_size, sig_dist, img_fname):

    print("Finding significant elements and plotting results ... "),
    sys.stdout.flush()

    # numpy.float64
    '''
    null_distn_params = {
        element:get_params(vals) for element,vals in null_distns.iteritems()
    }
    '''
    null_distn_params = {}
    to_del = []
    for element,vals in null_distns.iteritems():
        params = get_params(vals)
        if params["sd"] != 0:
            null_distn_params[element] = params
        else:
            to_del.append(element)

    for e in to_del:
        del null_distns[e]

    print(json.dumps(null_distn_params, sort_keys=True, indent=4))

    #float
    standard_null_distns = {
        element:standardise(vals, null_distn_params[element]) \
            for element,vals in null_distns.iteritems()
    }

    #float
    standard_query_values = {
        element:standardise([count], null_distn_params[element])[0] \
            for element,count in query_values.iteritems() if element in null_distn_params
    }

    #float
    query_distance = {
        element:abs(value) for element,value in standard_query_values.iteritems()
    }

    #float
    significant_elements = {
        element:distance for element,distance in query_distance.iteritems() \
            if distance >= sig_dist
    }

    draw_plot(
        standard_null_distns, standard_query_values, query_distance,
        null_distn_params, sample_size, sig_dist, img_fname
    )

    print("Done")
    sys.stdout.flush()

    return significant_elements

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("consensus_rti_intersect_fname")
    parser.add_argument("content_classified_rtis_fname")
    parser.add_argument("rmsk_align_bed_fname")
    parser.add_argument("sample_size", type=int)
    parser.add_argument("num_samples", type=int)
    parser.add_argument("min_prop", type=float)
    parser.add_argument("num_procs", type=int)
    parser.add_argument("consensus_tss_rmsk_align_intersect_fname")
    parser.add_argument("sig_dist", type=float)
    parser.add_argument("img_fname")
    parser.add_argument("element_class_json")
    args = parser.parse_args()

    with open(args.element_class_json, 'r') as f:
        element_class_dict = json.load(f)

    transcript_feat_distns \
        = get_transcript_feat_distns(args.consensus_rti_intersect_fname)

    rtis_by_classification \
        = get_rtis_by_classification(args.content_classified_rtis_fname)

    null_distns = generate_null_distns(
        transcript_feat_distns, rtis_by_classification,
        args.rmsk_align_bed_fname, args.sample_size, args.num_samples,
        args.min_prop, args.num_procs, element_class_dict
    )

    query_values = get_intersect_summary(
        args.consensus_tss_rmsk_align_intersect_fname, 10, args.min_prop, element_class_dict
    )
    print(query_values)

    significant_elements = plot_significant_elements(
        query_values, null_distns,
        args.sample_size, args.sig_dist, args.img_fname
    )

    # TODO: save results in a JSON file
