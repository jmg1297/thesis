'''
Given an mBED file of MSTRG and Ensembl transcripts, and a gffcompare tmap file,
find the StringTie transcripts in the tmap matching each of the Ensembl IDs.
'''

def get_tmap_dict(tmap_fname):
    tmap_dict = {}
    with open(tmap_fname, 'r') as f:
        headers = f.readline().strip().split()
        for line in f:
            line_list = line.strip().split()
            line_dict = {h:e for h,e in zip(headers, line_list)}
            if line_dict["class_code"] in ["=", "c", "j", "e"]:
                try:
                    tmap_dict[line_dict["ref_id"]]\
                        .append([line_dict["qry_id"], line_dict["class_code"]])
                except KeyError:
                    tmap_dict[line_dict["ref_id"]] \
                        = [[line_dict["qry_id"], line_dict["class_code"]]]
    return tmap_dict

def append_strg(mbed_fname, tmap_fname, output_fname):
    tmap_dict = get_tmap_dict(tmap_fname)

    with open(mbed_fname, 'r') as f, open(output_fname, 'wa') as out:
        for line in f:
            line_list = line.strip().split()
            ensembl_id = line_list[9]
            try:
                strg_matches = ";".join(
                    map(lambda x: ",".join(x), tmap_dict[ensembl_id]))
            except KeyError:
                strg_matches = "NONE"
            out_line = line.strip() + "\t" + strg_matches + "\n"
            out.write(out_line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mbed_fname")
    parser.add_argument("tmap_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    append_strg(args.mbed_fname, args.tmap_fname, args.output_fname)
