'''
Given a BED-like file with each line a reconstructed transcript, use a
gffcompare tmap file to find the class code for each transcript, and the
matching reference transcript and gene, if it exists. Also use the Ensembl GTF
file to add the biotype of matching reference transcripts.
'''

import sys
import re
import logging

logging.basicConfig(level=logging.INFO)

def get_tmap_dict(tmap_fname):
    logging.info("Parsing tmap file")
    tmap_dict = {}
    with open(tmap_fname, 'r') as f:
        headers = f.readline().strip().split()
        for line in f:
            line_dict = {h:v for h,v in zip(headers, line.strip().split())}
            tmap_dict[line_dict["qry_id"]] = {
                "ref_transcript_id":line_dict["ref_id"],
                "ref_gene_name":line_dict["ref_gene_id"],
                "class_code":line_dict["class_code"]}
    logging.info("Done")
    return tmap_dict

def get_attr_dict(attr_str):
    return dict(
        map(
            lambda y: re.split(' *"', y)[0:-1],
            map(
                lambda x: x.strip(),
                attr_str.split(";")[0:-1]
            )
        )
    )

def get_gtf_dict(gtf_fname):
    logging.info("Parsing GTF file")
    gtf_dict = {}
    with open(gtf_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split("\t")

            if line_list[2] == "transcript":
                attr_str = line_list[-1]
                attr_dict = get_attr_dict(attr_str)
                gtf_dict[attr_dict["transcript_id"]] = attr_dict["transcript_biotype"]
    logging.info("Done")
    return gtf_dict

def get_output_lines(bed_fname, tmap_dict, gtf_dict):
    logging.info("Parsing BED file")
    output_lines = []
    with open(bed_fname, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript_id = line_list[3]
            tmap_info = tmap_dict[transcript_id]
            biotype = gtf_dict.get(tmap_info["ref_transcript_id"], ".")

            output_line = "\t".join([
                line.strip(),
                tmap_info["class_code"],
                tmap_info["ref_transcript_id"],
                tmap_info["ref_gene_name"],
                biotype])
            output_lines.append(output_line)
    logging.info("Done")
    return output_lines

def write_output(output_lines, output_fname):
    logging.info("Writing output")
    with open(output_fname, 'wa') as out:
        out.write("\n".join(output_lines) + "\n")
    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("bed_fname")
    parser.add_argument("tmap_fname")
    parser.add_argument("gtf_fname")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    tmap_dict = get_tmap_dict(args.tmap_fname)

    gtf_dict = get_gtf_dict(args.gtf_fname)

    output_lines = get_output_lines(args.bed_fname, tmap_dict, gtf_dict)

    write_output(output_lines, args.output_fname)
