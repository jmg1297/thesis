'''
Given an mBED file of transcripts and RTIs, and a BED-like file of RTIs with
content classifications, write the each line of the mBED into a separate file
based on the RTI classification.
'''

def get_class_dict(rti_class_bed):
    class_dict = {}
    with open(rti_class_bed, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            class_dict[line_list[3]] = line_list[6]
    return class_dict

def split_mbed(mbed_fname, class_dict, output_base_fname):
    classes = set(class_dict.values())
    out_fname_dict = {c:".".join([output_base_fname, c, "mbed"]) for c in classes}
    output_dict = {c:[] for c in classes}

    with open(mbed_fname, 'r') as f:
        for line in f:
            rti_name = line.strip().split()[9]
            output_dict[class_dict[rti_name]].append(line)

    for c in classes:
        with open(out_fname_dict[c], 'wa') as f:
            f.write("".join(output_dict[c]))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("rti_class_bed")
    parser.add_argument("mbed_fname")
    parser.add_argument("output_base_fname")
    args = parser.parse_args()

    class_dict = get_class_dict(args.rti_class_bed)

    split_mbed(args.mbed_fname, class_dict, args.output_base_fname)
