'''
Run a Fisher's test on parents of retrocopies expressed in BL6 and conserved in
CAST
'''

from scipy.stats import fisher_exact
from numpy import array

def read_list(fname):
    with open(fname, 'r') as f:
        parent_list = [l.strip() for l in f]
    return parent_list

def get_counts(expressed_fname, conserved_fname, all_fname):
    '''
    Each input file should be a simple .txt list of parent transcript ENSEMBL IDs
    '''

    expressed = set(read_list(expressed_fname))
    conserved = set(read_list(conserved_fname))
    all_parents = set(read_list(all_fname))

    count_dict = {
        es:{
            cs:0 for cs in ["conserved", "!conserved"]
        } for es in ["expressed", "!expressed"]
    }

    count_dict["expressed"]["conserved"] = len(expressed.intersection(conserved))

    count_dict["expressed"]["!conserved"] = len(expressed.difference(conserved))

    count_dict["!expressed"]["conserved"] = len(conserved.difference(expressed))

    count_dict["!expressed"]["!conserved"] = len(all_parents) - sum(reduce(lambda a,d: a+d.values(), count_dict.values(), []))

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
            "Fisher's exact test on BL6 expressed retrocopy parents conserved in CAST\n"
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
    parser.add_argument("expressed_fname")
    parser.add_argument("conserved_fname")
    parser.add_argument("all_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    count_array = get_counts(args.expressed_fname, args.conserved_fname, args.all_fname)

    run_fisher(count_array, args.output_fname)
