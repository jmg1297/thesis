'''
Get the counts needed to run a Fisher's Exact Test to check whether retrocopy
expression in BL6 and conservation in CAST are related
'''

import subprocess
from scipy.stats import fisher_exact
from numpy import array

def get_counts(expressed_bed, conserved_bed, all_retrocopies_bed):
    count_dict = {
        es:{
            cs:0 for cs in ["conserved", "!conserved"]
        } for es in ["expressed", "!expressed"]
    }

    expressed = []
    with open(expressed_bed, 'r') as f:
        expressed = [l.strip().split()[3] for l in f]
    expressed = set(expressed)

    conserved = []
    with open(conserved_bed, 'r') as f:
        conserved = reduce(lambda a,x: a+x, [l.strip().split()[3].split(",") for l in f], [])
    conserved = set(conserved)

    all_rcs = []
    with open(all_retrocopies_bed, 'r') as f:
        all_rcs= [l.strip().split()[3] for l in f]
    all_rcs = set(all_rcs)

    count_dict["expressed"]["conserved"] = len(expressed.intersection(conserved))

    count_dict["expressed"]["!conserved"] = len(expressed.difference(conserved))

    count_dict["!expressed"]["conserved"] = len(conserved.difference(expressed))

    count_dict["!expressed"]["!conserved"] = len(all_rcs) - sum(reduce(lambda a,d: a+d.values(), count_dict.values(), []))

    count_array = array([
        [count_dict["expressed"]["conserved"], count_dict["expressed"]["!conserved"]],
        [count_dict["!expressed"]["conserved"], count_dict["!expressed"]["!conserved"]]
    ])

    return count_array

def run_fisher(count_array, output_fname):

    alts = ["two-sided", "less", "greater"]
    results_dict = {}
    for a in alts:
        results_dict[a] = fisher_exact(count_array, alternative=a)
    #odds_ratio, p_value = fisher_exact(count_array)

    with open(output_fname, 'wa') as f:
        f.write(
            "Fisher's exact test BL6 expressed retrocopies conserved in CAST"
        )
        f.write("Input matrix:\n")
        f.write(" "*14 + "conserved\n")
        f.write("\n")
        f.write(" "*15 + "+   -\n")
        f.write("\n")
        f.write(" "*12 + "+ {} {}\n".format(count_array[0][0], count_array[0][1]))
        f.write("expressed\n")
        f.write(" "*12 + "- {} {}\n".format(count_array[1][0], count_array[1][1]))
        f.write("\n")
        for a, (odds_ratio, p_value) in results_dict.items():
            f.write("Alternative: {}\n".format(a))
            f.write("\n")
            f.write("Odds ratio: {}\n".format(odds_ratio))
            f.write("p-value: {}\n".format(p_value))
            f.write("\n\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("expressed_bed")
    parser.add_argument("conserved_bed")
    parser.add_argument("all_retrocopies_bed")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    count_array = get_counts(args.expressed_bed, args.conserved_bed, args.all_retrocopies_bed)

    run_fisher(count_array, args.output_fname)
