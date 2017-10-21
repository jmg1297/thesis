'''
Split an extended mBED file based on whether retrotransposon content of
pairs of transcripts agrees
'''

import tempfile as tf
import subprocess as sub

def split_by_rt_agreement(mbed_fname, output_base_fname):
    output_dict = {"rt_agree":[], "rt_disagree":[], "no_rts":[]}
    out_fname_dict = {
        "rt_agree":output_base_fname + ".rt_agree.mbed",
        "rt_disagree":output_base_fname + ".rt_disagree.mbed",
        "no_rts":output_base_fname + ".no_rts.mbed"}

    with open(mbed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()

            cell_rt_info_str = line_list[6]
            if cell_rt_info_str != "NONE":
                cell_rt_info = set(cell_rt_info_str.split(","))
            else:
                cell_rt_info = "NONE"

            ensembl_strg_rt_info_str = line_list[13]
            if ensembl_strg_rt_info_str != "NONE":
                ensembl_strg_rt_info = set(
                    reduce(
                        lambda a,x: a+x,
                        [e.split(",")[2:]
                            for e in ensembl_strg_rt_info_str.split(";")],
                        []))
            else:
                ensembl_strg_rt_info = "NONE"

            if cell_rt_info == "NONE":
                if ensembl_strg_rt_info == "NONE":
                    output_dict["no_rts"].append(line)
                else:
                    output_dict["rt_disagree"].append(line)
            elif len(cell_rt_info.intersection(ensembl_strg_rt_info)) == 0:
                output_dict["rt_disagree"].append(line)
            else:
                output_dict["rt_agree"].append(line)

    for k,l in output_dict.items():
        with tf.NamedTemporaryFile() as tmp:
            tmp.write("".join(l))
            tmp.flush()
            sub.call(
                "sort -k1,1 -k2,2n {} > {}".format(tmp.name, out_fname_dict[k]),
                shell=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mbed_fname")
    parser.add_argument("output_base_fname")
    args = parser.parse_args()

    split_by_rt_agreement(args.mbed_fname, args.output_base_fname)
