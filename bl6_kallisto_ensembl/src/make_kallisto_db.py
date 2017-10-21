'''
Import kallisto quant results into an SQLite database and create tables with
cell-type/sex expression averages for both transcripts and genes
'''

import sys
import os
import re
import sqlite3
from numpy import mean

def scrub(s):
    return re.sub("[\";]", "", s)

def get_gene_transcript_dict(annotation_fname):
    print("Getting transcript to gene dictionary ... "),
    sys.stdout.flush()

    gene_transcript_dict = {}
    with open(annotation_fname, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue

            line_list = line.strip().split()
            if line_list[2] != "transcript":
                continue

            gene_id = scrub(line_list[9])
            gene_name = scrub(line_list[17])
            transcript_id = scrub(line_list[13])
            gene_transcript_dict[transcript_id] \
                = {"gene_id":gene_id, "gene_name":gene_name}
    print("Done")
    return gene_transcript_dict

def make_samples_table(db_fname, sample_dir, quant_fname, gene_transcript_dict):
    print("Making samples table ... ")
    sys.stdout.flush()

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    transcript_data = {}
    gene_data = {}
    tdata = []

    for sample in os.listdir(sample_dir):
        print("Making table for sample {}".format(sample))
        sys.stdout.flush()

        code = "_".join(sample.split("_")[0:2])

        filepath = os.path.join(sample_dir, sample, quant_fname)
        with open(filepath, 'r') as f:
            headers = f.readline().strip().split()
            for line in f:
                line_list = line.strip().split()
                transcript_id = line_list[0]
                length = int(line_list[1])
                eff_length, est_counts, tpm = map(float, line_list[2:])
                tdata.append((code, transcript_id, length, eff_length, est_counts, tpm))
                try:
                    transcript_data[transcript_id][code] = tpm
                except KeyError:
                    transcript_data[transcript_id] = {code:tpm}

                try:
                    gene_info = gene_transcript_dict[transcript_id]
                except KeyError:
                    tmp_id = transcript_id.split("_")[0]
                    try:
                        gene_info = gene_transcript_dict[tmp_id]
                    except KeyError:
                        print("No matching gene for transcript {}".format(transcript_id))
                        continue
                gene_id = gene_info['gene_id']
                gene_name = gene_info['gene_name']
                try:
                    gene_data[gene_info["gene_id"]]["tpms"][code] += tpm
                except KeyError:
                    if gene_id not in gene_data:
                        gene_data[gene_info["gene_id"]] \
                            = {"gene_name":gene_info["gene_name"], "tpms":{code:tpm}}
                    else:
                        gene_data[gene_id]["tpms"][code] = tpm

    cur.execute(
        '''
        CREATE TABLE IF NOT EXISTS transcripts (
        sample TEXT NOT NULL,
        transcript_id TEXT NOT NULL,
        length INT NOT NULL,
        eff_length REAL NOT NULL,
        est_counts REAL NOT NULL,
        tpm REAL NOT NULL
        )
        '''
    )
    seq = ",".join(["?"]*len(tdata[0]))
    cur.executemany("INSERT INTO transcripts VALUES ({})".format(seq), tdata)

    cur.execute(
        '''
        CREATE TABLE IF NOT EXISTS genes (
        sample TEXT NOT NULL,
        gene_id TEXT NOT NULL,
        gene_name TEXT NOT NULL,
        tpm REAL NOT NULL
        )
        '''
    )

    gdata = []
    for gene_id, gene_info in gene_data.iteritems():
        gene_name = gene_info["gene_name"]
        tpms = gene_info["tpms"]
        for code,tpm in tpms.iteritems():
            gdata.append((code, gene_id, gene_name, tpm))

    cur.executemany(
        "INSERT INTO genes VALUES ({})".format(",".join(["?"]*len(gdata[0]))),
        gdata
    )

    conn.commit()

    print("Finished samples table")
    return transcript_data, gene_data

def make_avgs_tables(db_fname, transcript_data, gene_data, metadata_db):
    '''
    Make one table containing gene expression values, averaged across cell-type,
    sex, sex+cell, all. Same again for transcripts.
    '''
    print("Making averages tables ... "),
    sys.stdout.flush()

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()

    mdconn = sqlite3.connect(metadata_db)
    mdcur = mdconn.cursor()

    query = '''
    SELECT DISTINCT(cell) FROM metadata WHERE hybrid=0 AND strain='C57BL-6J'
    '''
    cell_types = [r[0] for r in mdcur.execute(query).fetchall()]

    query = '''
    SELECT DISTINCT(sex) FROM metadata WHERE hybrid=0 AND strain='C57BL-6J'
    '''
    sexes = [r[0] for r in mdcur.execute(query).fetchall()]

    cell_samples = {}
    for cell in cell_types:
        query = '''
        SELECT code FROM metadata
        WHERE hybrid=0 AND strain='C57BL-6J' AND cell=?
        '''
        cell_samples[cell] = [r[0] for r in mdcur.execute(query, (cell,)).fetchall()]

    sex_samples = {}
    for sex in sexes:
        query = '''
        SELECT code FROM metadata
        WHERE hybrid=0 AND strain='C57BL-6J' AND sex=?
        '''
        sex_samples[sex] = [r[0] for r in mdcur.execute(query, (sex,)).fetchall()]

    cross_samples = {}
    for sex in sexes:
        for cell in cell_types:
            query = '''
            SELECT code FROM metadata
            WHERE hybrid=0 AND strain='C57BL-6J' AND sex=? AND cell=?
            '''
            cross_samples[(sex,cell)] = [r[0] for r in mdcur.execute(query, (sex, cell)).fetchall()]

    query = '''
    SELECT code FROM metadata
    WHERE hybrid=0 AND strain='C57BL-6J'
    '''
    all_samples = [r[0] for r in mdcur.execute(query).fetchall()]

    cur.execute(
        '''
        CREATE TABLE IF NOT EXISTS transcript_averages (
        transcript_id TEXT,
        b_cell REAL,
        t_cell REAL,
        female REAL,
        male REAL,
        female_b REAL,
        female_t REAL,
        male_b REAL,
        male_t REAL,
        all_samples REAL
        )
        '''
    )

    '''
    transcript_id TEXT,
    b_cell REAL,
    t_cell REAL,
    female REAL,
    male REAL,
    female_b REAL,
    female_t REAL,
    male_b REAL,
    male_t REAL,
    all REAL
    '''

    data = []
    for transcript_id, tpm_dict in transcript_data.iteritems():
        b_avg = mean([tpm_dict[sample] for sample in cell_samples['B']])
        t_avg = mean([tpm_dict[sample] for sample in cell_samples['T']])
        f_avg = mean([tpm_dict[sample] for sample in sex_samples['Female']])
        m_avg = mean([tpm_dict[sample] for sample in sex_samples['Male']])
        fb_avg = mean([tpm_dict[sample] for sample in cross_samples[('Female', 'B')]])
        ft_avg = mean([tpm_dict[sample] for sample in cross_samples[('Female', 'T')]])
        mb_avg = mean([tpm_dict[sample] for sample in cross_samples[('Male', 'B')]])
        mt_avg = mean([tpm_dict[sample] for sample in cross_samples[('Male', 'T')]])
        all_avg = mean([tpm_dict[sample] for sample in all_samples])
        data.append(
            (transcript_id, b_avg, t_avg, f_avg, m_avg, fb_avg, ft_avg, mb_avg, mt_avg, all_avg)
        )
    cur.executemany(
        '''
        INSERT INTO transcript_averages VALUES ({})
        '''.format(",".join(["?"]*len(data[0]))),
        data
    )

    conn.commit()

    cur.execute(
        '''
        CREATE TABLE IF NOT EXISTS gene_averages (
        gene_id TEXT,
        gene_name TEXT,
        b_cell REAL,
        t_cell REAL,
        female REAL,
        male REAL,
        female_b REAL,
        female_t REAL,
        male_b REAL,
        male_t REAL,
        all_samples REAL
        )
        '''
    )

    data = []
    for gene_id, gene_info in gene_data.iteritems():
        tpm_dict = gene_info["tpms"]
        print(tpm_dict.keys())
        gene_name = gene_info["gene_name"]
        b_avg = mean([tpm_dict[sample] for sample in cell_samples['B']])
        t_avg = mean([tpm_dict[sample] for sample in cell_samples['T']])
        f_avg = mean([tpm_dict[sample] for sample in sex_samples['Female']])
        m_avg = mean([tpm_dict[sample] for sample in sex_samples['Male']])
        fb_avg = mean([tpm_dict[sample] for sample in cross_samples[('Female', 'B')]])
        ft_avg = mean([tpm_dict[sample] for sample in cross_samples[('Female', 'T')]])
        mb_avg = mean([tpm_dict[sample] for sample in cross_samples[('Male', 'B')]])
        mt_avg = mean([tpm_dict[sample] for sample in cross_samples[('Male', 'T')]])
        all_avg = mean([tpm_dict[sample] for sample in all_samples])
        data.append(
            (gene_id, gene_name, b_avg, t_avg, f_avg, m_avg, fb_avg, ft_avg, mb_avg, mt_avg, all_avg)
        )
    cur.executemany(
        '''
        INSERT INTO gene_averages VALUES ({})
        '''.format(",".join(["?"]*len(data[0]))),
        data
    )

    conn.commit()
    print("Done!")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation_fname")
    parser.add_argument("metadata_db")
    parser.add_argument("db_fname")
    parser.add_argument("sample_dir")
    parser.add_argument("quant_fname")
    args = parser.parse_args()

    gene_transcript_dict = get_gene_transcript_dict(args.annotation_fname)

    transcript_data, gene_data = make_samples_table(
        args.db_fname, args.sample_dir, args.quant_fname, gene_transcript_dict
    )

    make_avgs_tables(
        args.db_fname, transcript_data, gene_data, args.metadata_db
    )
