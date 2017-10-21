'''
Given a CSV file of StringTie merge transcripts and Ensembl transcript IDs,
and TPM data for StringTie merge and Ensembl transcripts for a set of samples,
draw a scatter plot of StringTie merge TPMs vs Ensembl TPMs.
'''

import os
import sys
import logging
import traceback as tb

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn import linear_model
from scipy.stats import pearsonr, spearmanr

logging.basicConfig(level=logging.INFO)

def parse_tpm_file(tpm_fname, log_tpm):
    logging.debug("Parsing TPM file {}".format(tpm_fname))
    tpms = {}
    with open(tpm_fname, 'r') as f:
        headers = f.readline().strip().split()
        for line in f:
            line_list = line.strip().split()
            transcript = line_list[0]
            if log_tpm:
                tpm = np.log10(float(line_list[4])+1)
            else:
                tpm = float(line_list[4])
            tpms[transcript] = tpm
    logging.debug("Finished")
    return tpms

def get_tpm_dict(ensembl_tpm_dir, strg_merge_tpm_dir, log_tpm, minimum_tpm):

    logging.info("Getting TPM values")

    ensembl_samples = os.listdir(ensembl_tpm_dir)
    strg_samples = os.listdir(strg_merge_tpm_dir)

    sample_list = list(set(ensembl_samples).intersection(set(strg_samples)))

    logging.info("{} samples found in both directories".format(len(sample_list)))

    ensembl_tpm_fnames = [os.path.join(ensembl_tpm_dir, s, "abundance.tsv")
        for s in sample_list]
    strg_tpm_fnames = [os.path.join(strg_merge_tpm_dir, s, "abundance.tsv")
        for s in sample_list]

    tpm_dict = {
        s:{"ensembl":parse_tpm_file(e, log_tpm), "strg":parse_tpm_file(m, log_tpm)}
        for s,e,m in zip(sample_list, ensembl_tpm_fnames, strg_tpm_fnames)}

    logging.info("Finished getting TPM values")
    return tpm_dict, sample_list

def get_plot_data(csv_fname, sample_list, tpm_dict, log_tpm, minimum_tpm):
    logging.info("Getting plot data")
    # Expected order is StringTie,Ensembl
    transcript_pairs = [line.strip().split(",") for line in open(csv_fname, 'r')]

    logging.info(
        "{} transcript pairs found, getting TPM values"
            .format(len(transcript_pairs)))

    plot_data = {s:{"strg":[], "ensembl":[]} for s in sample_list}

    if log_tpm:
        minimum_tpm = np.log10(minimum_tpm + 1)

    for s in sample_list:
        logging.debug("Processing sample {}".format(s))
        for strg,ensembl in transcript_pairs:
            try:
                stpm = tpm_dict[s]["strg"][strg]
                etpm = tpm_dict[s]["ensembl"][ensembl]
                if stpm > minimum_tpm and etpm > minimum_tpm:
                    plot_data[s]["strg"].append(stpm)
                    plot_data[s]["ensembl"].append(etpm)
            except KeyError:
                # Probably a transcript on a scaffold, therefore not in the
                # fasta file
                continue

    logging.info("Finished getting plot data")
    return plot_data

def get_color_dict(sample_list):
    # Really specific code, unfortunately - apologies future users
    cmap = plt.get_cmap('tab20c')
    colors = [cmap(i/20.0) for i in range(16) if (i+1)%4 != 0]
    color_dict = {s:c for s,c in zip(sorted(sample_list), colors)}
    return color_dict

'''
def run_linear_regr(plot_data, sample_list):
    logging.info("Running linear model")
    #X = np.array([plot_data[s]["strg"] for s in sorted(sample_list)])
    #y = np.array([plot_data[s]["ensembl"] for s in sorted(sample_list)])

    X = np.array(
        reduce(
            lambda a,x: a+x,
            [plot_data[s]["strg"] for s in sorted(sample_list)[0:1]],
            []))

    y = np.array(
        reduce(
            lambda a,x: a+x,
            [plot_data[s]["ensembl"] for s in sorted(sample_list)[0:1]],
            []))

    a = np.array([plot_data[s]["strg"] for s in sorted(sample_list)[0:1]]).T
    b = np.array([plot_data[s]["ensembl"] for s in sorted(sample_list)[0:1]]).T

    #regr = linear_model.LinearRegression()
    #regr.fit(X, y)

    #pred = regr.predict(X)
    #logging.info("Finished")
    #return regr.coef_, pred

    pearson_coeff, pearson_p = pearsonr(X, y)

    print("Pearson coeff: {}".format(pearson_coeff))
    print("Pearson p: {}".format(pearson_p))

    spearman_coeff, spearman_p = spearmanr(X, y)

    print("Spearman coeff: {}".format(spearman_coeff))
    print("Spearman p: {}".format(spearman_p))
'''

def get_corr_coef(x, y):
    coeff, p = spearmanr(x, y)
    return coeff, p

def draw_scatter(plot_data, sample_list, log_tpm, img_fname):
    logging.info("Drawing scatter plot")
    color_dict = get_color_dict(sample_list)

    fig,ax = plt.subplots(1,1)
    fig.set_size_inches(10,10)

    rho = []
    p_vals = []

    for sample in reversed(sorted(sample_list)):
        tpm_dict = plot_data[sample]

        x = tpm_dict["strg"]
        y = tpm_dict["ensembl"]

        color = color_dict[sample]
        ax.scatter(x, y, c=color, s=6, alpha=0.6, lw=0, label=sample)

        corr_coef, p = get_corr_coef(x,y)
        rho.append(corr_coef)
        p_vals.append(p)


        ax.plot(
            np.unique(x),
            np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),
            color=color,
            linestyle='dashed',
            lw=2)
    ax.set_title("Spearman's rho = %.2f, p = %.2f" % (np.mean(rho), np.mean(p_vals)))

    #print("Linear model coefficients: {}".format(lm_coef))
    #print(lm_pred)

    lim = max(ax.get_xlim() + ax.get_ylim())
    ax.set_xlim([-0.05,lim])
    ax.set_ylim([-0.05,lim])

    if log_tpm:
        xlabel = "StringTie $\log_{10}(TPM+1)$"
        ylabel = "Ensembl $\log_{10}(TPM+1)$"
    else:
        xlabel = "StringTie TPM"
        ylabel = "Ensembl TPM"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.legend()

    logging.info("Finished scatter plot, saving figure to {}".format(img_fname))
    plt.savefig(img_fname)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_fname")
    parser.add_argument("ensembl_tpm_dir")
    parser.add_argument("strg_merge_tpm_dir")
    parser.add_argument("-l", "--log_tpm", type=int, default=1)
    parser.add_argument("-m", "--minimum_tpm", type=float, default=0)
    parser.add_argument("img_fname")
    args = parser.parse_args()

    tpm_dict, sample_list = get_tpm_dict(
        args.ensembl_tpm_dir, args.strg_merge_tpm_dir,
        args.log_tpm, args.minimum_tpm)

    plot_data = get_plot_data(
        args.csv_fname, sample_list, tpm_dict, args.log_tpm, args.minimum_tpm)

    draw_scatter(
        plot_data, sample_list, args.log_tpm, args.img_fname)
