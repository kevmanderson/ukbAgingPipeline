#!/bin/python

import os
import json
import shutil
import subprocess
import argparse
import sqlite3
import glob
import pandas as pd

def create_out_dir_name(imp_or_cal, maf, hwe, mind, geno, info):
    '''
    Stitch together preprocessing parameters into a descriptive folder name
    :param maf: Minor-Allele Frequency, e.g. 0.05
    :param hwe: Hardy-Weinberg Equilibrium, e.g. 1e-6
    :param mind: SNP Missingness
    :param geno: Genotype Missingness
    :param info: Imputation Quality Confidence Score
    :return:
    '''
    out_dir = 'ukb_{}_v3_maf{}_hwe{}_mind{}_geno{}_info{}'.format(imp_or_cal, maf, hwe, mind, geno, info)
    out_dir = out_dir.replace('.', '_').replace('e-', 'e')
    return out_dir

def create_subject_filter_list(genetic_sex_df, sex_df, aneuploidy_df, zygosity_df, ethnicity_df, kinship_df):
    # start with subjects that have a genetic sex report
    keep_subs = list(genetic_sex_df['eid'])

    # :::: SUBJECT FILTER ::::
    # GENETIC vs REPORTED SEX
    genetic_sex_df['genetic_sex'] = ['Male' if float(x) == 1 else 'Female' for x in genetic_sex_df.value]

    # genetic ethnicity
    sex_df['sex'] = ['Male' if float(x) == 1 else 'Female' for x in sex_df.value]
    sex_df = sex_df.merge(genetic_sex_df, on='eid')

    # only keep subjects with reported sex matching genetic sex
    mismatch_sex_df = sex_df.loc[sex_df['sex'] != sex_df['genetic_sex']]
    keep_subs = set(keep_subs) - set(mismatch_sex_df['eid'])

    # :::: SUBJECT FILTER ::::
    # SEX ANEUPLOIDY
    keep_subs = set(keep_subs) - set(aneuploidy_df['eid'])

    # :::: SUBJECT FILTER ::::
    # European ancestry
    keep_subs = keep_subs.intersection(ethnicity_df['eid'])

    # :::: SUBJECT FILTER ::::
    # Excluded from kinship
    rm_kinship_df = kinship_df.loc[kinship_df['value'].astype(float) == -1]
    keep_subs = keep_subs - set(rm_kinship_df['eid'])

    return list(keep_subs)


def preprocess_bgen_files(chrom, maf, hwe, mind, geno, sub_keep_file, imputed_write_dir, plink2):
    '''

    Parameters
    ----------
    pheno: str, required
        The UKB phenotype to be predicted. e.g. fluidint, height, neuroticism
    chrom: str, required
        1-22
    maf: str, required
        Minor Allele Frequency (MAF). SNPs with MAF lower than this threshold will be excluded.
    hwe: str, required
        Hardy-Weinberg Equilibrium (HWE). p-value cutoff to eliminate SNPs not in HWE.
    mind: str, required
        Missingness (subject-wise). Remove subjects with SNP missingness above this value
    geno: str, required
        Genotype Rate. Censor SNPs with genotype rates lower than this value
    scratch_dir: str, required
        Root folder where the big data files are being stores


    Returns
    ----------
    string
        returns the 'plink_submit' object, corresponding to the PLINK command to be executed.

        This function will write and submit the PLINK preprocessing command to our SLURM server

    '''

    # path to processed bgen file
    cur_bgen    = '/gpfs/milgram/data/UKB/ukb_snp/ukb_imp_chr{}_v3.bgen'.format(chrom)
    cur_sample  = '/gpfs/milgram/data/UKB/ukb_snp/ukb25163_imp_chr{}_v3_s487324.sample'.format(chrom)

    # new outputs
    bgen_out    = os.path.join(imputed_write_dir, 'ukb_eur_imp_chr{}_v3_info090_maf{}_hwe{}_mind{}_geno{}'.format(chrom, maf, hwe, mind, geno))
    bgen_out    = bgen_out.replace('0.0', '0').replace('.', '')

    # sample
    snp_file    = '/gpfs/milgram/scratch/holmes/buckner_aging/data/ukb/snp_info/ukb_mfi_chr{}_info90_plink_positions.txt'.format(chrom)

    # plink preprocessing command
    plink_cmd   = (f'''{plink2} --bgen {cur_bgen} ref-first --sample {cur_sample} --max-alleles 2 --maf {maf} --hwe {hwe} --mind {mind} --geno {geno} -export bgen-1.2 'bits=8' --out {bgen_out} --keep {sub_keep_file} --extract 'range' {snp_file}''') # --memory 120000

    # convert bgen to bed/bim/fam
    plink_cmd_2 = (f'''{plink2} --bgen {bgen_out}.bgen ref-first --sample {bgen_out}.sample --make-bed --out {bgen_out}''')

    # submit both at once
    plink_submit = plink_cmd + '\n\n' + plink_cmd_2
    plink_submit = plink_cmd_2

    slurm_out_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging/slurm/'
    slurm_file  = os.path.join(slurm_out_dir, 'plink_imp_chr{}'.format(chrom))
    cmd_path    = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=20, cmd=plink_submit, stime='06:00:00', jobName='imp_'+str(chrom))
    job_id      = submitSlurm(cmd_path, dependencies=None)

    return slurm_file

def retrieve_field(conn, field, table_name, visit):
    cur = conn.cursor()
    cur.execute("SELECT * FROM {} WHERE field=='{}'".format(table_name, field))
    colnames = [x[0] for x in cur.description]
    rows   = cur.fetchall()
    out_df = pd.DataFrame(rows)
    out_df.columns = colnames
    out_df = out_df.loc[out_df.time == visit]
    return out_df

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn

def main():
    example_text = '''example: 
     python test.py -t template/test.py'''

    parser = argparse.ArgumentParser(epilog=example_text, add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', '-c', dest='config', required=True, help='Full path to the user configuration file')
    parser.add_argument('--chrom', dest='maf', default=0.01, required=False, help='Minor Allele Frequency')
    parser.add_argument('--maf', dest='maf', default=0.01, required=False, help='Minor Allele Frequency')
    parser.add_argument('--hwe', dest='hwe', default=1e-6, required=False, help='Hardy-Weinberg Equilibrium')
    parser.add_argument('--mind', dest='mind', default=0.02, required=False, help='Genotype Missingness')
    parser.add_argument('--geno', dest='geno', default=0.02, required=False, help='SNP Missingness')
    parser.add_argument('--info', dest='info', default=0.60, required=False, help='SNP Missingness')
    parser.add_argument('--european', dest='info', default=0.60, required=False, help='SNP Missingness')
    parser.add_argument('--help', '-h', action='help', default=argparse.SUPPRESS, help='This script will unpack and convert UK Biobank')
    opt = parser.parse_args()
    config_file = opt.config

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.maf   = 0.01
        opt.hwe   = 1e-6
        opt.mind  = 0.02
        opt.geno  = 0.02
        opt.info  = 0.60
        opt.chrom = 22
        opt.slurm = True

    # read user configuration file
    with open(opt.config, 'r') as f:
        config_json = json.load(f)[0]

    # sqlite
    sqlite_path = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukbcc_sqlite3_new.db'
    conn   = create_connection(sqlite_path)

    # define directories
    imputed_dir = os.path.join(config_json['base_dir'], 'data/ukb/genotypes/imputed')
    imputed_out_folder  = create_out_dir_name(imp_or_cal='imp', maf=opt.maf, hwe=opt.hwe, mind=opt.mind, geno=opt.geno, info=opt.info)
    called_out_folder   = create_out_dir_name(imp_or_cal='cal', maf=opt.maf, hwe=opt.hwe, mind=opt.mind, geno=opt.geno, info=opt.info)

    # create output directory
    proc_dir = os.path.join(config_json['base_dir'], 'data/ukb/genotypes/processed')
    if not os.path.exists(proc_dir):
        os.mkdir(proc_dir)
    imputed_out_path = os.path.join(proc_dir, imputed_out_folder)
    if not os.path.exists(imputed_out_path):
        os.mkdir(imputed_out_path)
    called_out_path = os.path.join(proc_dir, called_out_folder)
    if not os.path.exists(imputed_out_path):
        os.mkdir(imputed_out_path)

    # identify subjects to keep for genetic QC
    genetic_sex_df = retrieve_field(conn=conn, field='22001', table_name='str', visit=0)
    sex_df         = retrieve_field(conn=conn, field='31', table_name='str', visit=0)
    aneuploidy_df  = retrieve_field(conn=conn, field='22019', table_name='str', visit=0)
    ethnicity_df   = retrieve_field(conn=conn, field='22006', table_name='str', visit=0)
    kinship_df     = retrieve_field(conn=conn, field='22021', table_name='str', visit=0)
    keep_subjects  = create_subject_filter_list(genetic_sex_df, sex_df, aneuploidy_df, zygosity_df, ethnicity_df, kinship_df)

    # identify bi-allelic SNPs with
    print('Reading UKB MFI')
    mfi_df = pd.read_table(os.path.join(imputed_dir, 'ukb_mfi_chr{}_v3.txt').format(opt.chrom), header=None)

    # filter based on INFO score
    mfi_filter_df = mfi_df.loc[mfi_df[7] >= opt.info]

    # only retain bi-allelic SNPs
    mfi_filter_df['a1_len'] = [len(x) for x in mfi_filter_df[3]]
    mfi_filter_df['a2_len'] = [len(x) for x in mfi_filter_df[4]]
    mfi_filter_df = mfi_filter_df.loc[mfi_filter_df['a1_len'] == 1]
    mfi_filter_df = mfi_filter_df.loc[mfi_filter_df['a2_len'] == 1]
    mfi_filter_df['chr'] = opt.chrom

    # write PLINK SNP filter file
    snp_filter_file = os.path.join(imputed_dir, 'ukb_plink_snp_filter_{}_v3.txt'.format(opt.chrom))
    mfi_filter_df[['chr', 2, 2]].to_csv(snp_filter_file, sep='\t', index=None, header=None)

    # write PLINK SNP subject filter
    subj_filter_file = os.path.join(imputed_dir, 'ukb_plink_subject_filter.txt')
    plink_df = pd.DataFrame({'FID': list(keep_subs), 'IID': list(keep_subs)})
    plink_df.to_csv(subj_filter_file, sep='\t', index=None)


    # plink preprocessing command
    plink2     = '/opt/plink2/plink2'
    cur_bgen   = os.path.join(imputed_dir, 'ukb_imp_chr{}_v3.bgen'.format(opt.chrom))
    cur_sample = os.path.join(imputed_dir, 'ukb25163_imp_chr{}_v3_s487324.sample'.format(opt.chrom))
    bgen_out   = os.path.join(imputed_dir, 'TMP_filter_ukb_imp_chr{}_v3.bgen'.format(opt.chrom))

    plink_cmd  = (f'''{plink2} \\\n--bgen {cur_bgen} ref-first \\\n--sample {cur_sample} \\\n--max-alleles 2 \\\n-export bgen-1.2 'bits=8' \\\n--out {bgen_out} \\\n--keep {subj_filter_file} \\\n--extract 'range' {snp_filter_file}''') # --memory 120000
    print(plink_cmd)

    mfi_list = [pd.read_table(x, header=None) for x in glob.glob(os.path.join(imputed_dir, 'ukb_mfi_chr*_v3.txt'))]
    mfi_df   = pd.concat(mfi_list)
    keep_subjects = create_subject_filter_list(genetic_sex_df, sex_df, aneuploidy_df, zygosity_df, ethnicity_df, kinship_df)


    # genotyped

    # subset to European

    # SNP QC

    # Genetic Relatedness





if __name__ == "__main__":
    main()









