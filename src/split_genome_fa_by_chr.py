import os
import sys
import argparse

def split_fasta(input_fasta, output_dir):
    output_file = None
    with open(input_fasta, 'r') as f:
        for line in f:
            if line[0] == ">":
                try:
                    output_file.close()
                except AttributeError:
                    pass
                print(line.strip())
                sys.stdout.flush()
                output_fname = os.path.join(output_dir, line.strip()[1:]+".fa")
                output_file = open(output_fname, "wa")
                output_file.write(line)
            else:
                output_file.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    split_fasta(args.input_fasta, args.output_dir)
