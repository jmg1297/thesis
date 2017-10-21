'''
Script to carry out PCA on a set of alignments and produce a plot showing the
first two principal components. Also do K-means clustering and use clusters to
color points on the plot.
'''

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn import preprocessing

def get_data(aln_out_fname, query_feat_fname):

    aln_dict = {}
    with open(aln_out_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            try:
                aln_dict[line_list[0]].append([
                    line_list[1],
                    float(line_list[2]),
                    int(line_list[3]),
                    int(line_list[5]),
                    float(line_list[10]),
                    float(line_list[11])
                ])
            except KeyError:
                aln_dict[line_list[0]] = [[
                    line_list[1],
                    float(line_list[2]),
                    int(line_list[3]),
                    int(line_list[5]),
                    float(line_list[10]),
                    float(line_list[11])
                ]]

    feat_length_dict = {}
    with open(query_feat_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            feat_length_dict[line_list[3]] = {
                "chrom":line_list[0],
                "length":int(line_list[2]) - int(line_list[1])
            }

    data = []

    f = open("high_aln.txt", 'wa')

    for rc, l in aln_dict.iteritems():
        for dl in l:
            if feat_length_dict[rc]["chrom"] == dl[0]:
                dl[2] = dl[2]/feat_length_dict[rc]["length"]
                p = ((dl[3])/100.0)*dl[2]
                dl.append(p)

                if dl[2] > 0.8:
                    f.write(rc + "\n")

                    data.append(dl[1:])

    f.close()
    return data

def plot_hists(data):
    fig, ax = plt.subplots(6, 1)
    fig.set_size_inches(4, 30)

    identity = []
    prop_len = []
    gap_opens = []
    evalue = []
    bit_score = []
    prop_bp = []

    for r in data:
        identity.append(r[0])
        prop_len.append(r[1])
        gap_opens.append(r[2])
        evalue.append(r[3])
        bit_score.append(r[4])
        prop_bp.append(r[5])

    #plt.yscale('log', nonposy='clip')

    for i,d in enumerate([identity, prop_len, gap_opens, evalue, bit_score, prop_bp]):
        ax[i].hist(d, bins=50)
        #ax[i].set_yscale('log', nonposy='clip')

    plt.savefig("aln_hists.svg")

def run_pca(data, out_summary_fname):
    pca = PCA()
    X = np.array(data)
    #X = X[:, [0,1]]
    transformed = pca.fit_transform(X)
    print(pca.components_)
    print(pca.explained_variance_ratio_)
    return transformed

def run_kmeans(data, out_summary_fname):
    X = np.array(data)
    #X = X[:, [0,1]]
    kmeans = KMeans(n_clusters=2)
    kmeans.fit(data)
    return kmeans.labels_

def plot_pca(pca_results, cluster_labels, plot_title, img_fname):
    plot_data = pca_results[:,:2]
    color_dict = {0:'red', 1:'blue'}
    colors = [color_dict[l] for l in cluster_labels]

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    ax.scatter(plot_data[:, 0], plot_data[:, 1], s=2, c=colors, lw=0, alpha=0.1)

    plt.savefig(img_fname)

def run_tsne(data):
    X = np.array(data)
    Xs = preprocessing.scale(X)
    model = TSNE(n_components=2, random_state=0)
    Y = model.fit_transform(X)
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.scatter(Y[:,0], Y[:,1], lw=0, c='black', alpha=0.2, s=4)
    plt.savefig("tsne.svg")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("aln_out_fname")
    parser.add_argument("query_feat_fname")
    parser.add_argument("out_summary_fname")
    parser.add_argument("plot_title")
    parser.add_argument("img_fname")
    args = parser.parse_args()

    data = get_data(args.aln_out_fname, args.query_feat_fname)
    #print(data)

    plot_hists(data)

    pca_results = run_pca(data, args.out_summary_fname)

    cluster_labels = run_kmeans(data, args.out_summary_fname)

    plot_pca(pca_results, cluster_labels, args.plot_title, args.img_fname)

    run_tsne(data)
