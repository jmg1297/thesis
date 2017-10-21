'''
Given: an extended mBED file with cell-specific StringTie transcripts,
associated Ensembl transcript IDs, and StringTie transcripts matching the
Ensembl transcripts; an mBED file of cell-specific StringTie transcripts with
overlapping retrotransposons; and an mBED file of ALL StringTie transcripts
with overlapping retrotransposon
then write out a file with retrotransposon information appended
'''

def get_rt_dict(transcript_rt_mbed):
    rt_dict = {}
    with open(transcript_rt_mbed, 'r') as f:
        for line in f:
            line_list = line.strip().split()
            transcript = line_list[3]
            rt = ":".join(line_list[9].split(":")[0:2])
            try:
                rt_dict[transcript].append(rt)
            except KeyError:
                rt_dict[transcript] = [rt]
    return rt_dict

def append_rt_overlaps(mbed_fname, cell_strg_rt_mbed, all_strg_rt_mbed, output_fname):
    cell_rt_dict = get_rt_dict(cell_strg_rt_mbed)
    all_rt_dict = get_rt_dict(all_strg_rt_mbed)

    with open(mbed_fname, 'r') as f, open(output_fname, 'wa') as out:
        for line in f:
            line_list = line.strip().split()
            cell_transcript = line_list[3]
            if line_list[-1] != "NONE":
                all_transcripts = [e.split(",") for e in line_list[-1].split(";")]
            else:
                all_transcripts = "NONE"

            try:
                cell_rt_info = ",".join(list(set(cell_rt_dict[cell_transcript])))
            except KeyError:
                cell_rt_info = "NONE"

            all_transcript_info = []
            if all_transcripts == "NONE":
                all_transcript_info_str = "NONE"
            else:
                for transcript,class_code in all_transcripts:
                    try:
                        rt_info = ",".join(list(set(all_rt_dict[transcript])))
                    except KeyError:
                        rt_info = "NONE"

                    info = ",".join([transcript, class_code, rt_info])
                    all_transcript_info.append(info)
                all_transcript_info_str = ";".join(all_transcript_info)

            output_line = "\t".join(
                line_list[0:6] \
                + [cell_rt_info] \
                + line_list[6:-1] \
                + [all_transcript_info_str])

            out.write(output_line + "\n")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mbed_fname")
    parser.add_argument("cell_strg_rt_mbed")
    parser.add_argument("all_strg_rt_mbed")
    parser.add_argument("output_fname")
    args = parser.parse_args()

    append_rt_overlaps(
        args.mbed_fname,
        args.cell_strg_rt_mbed, args.all_strg_rt_mbed,
        args.output_fname)
