#!/bin/python

import os
import json
import shutil
import subprocess
import argparse
import sqlite3
import glob
import pandas as pd

def snp_dirname(imp_or_cal, maf, hwe, mind, geno, info):
    '''
    Stitch together preprocessing parameters into a descriptive folder name
    :param maf: Minor-Allele Frequency, e.g. 0.05
    :param hwe: Hardy-Weinberg Equilibrium, e.g. 1e-6
    :param mind: SNP Missingness
    :param geno: Genotype Missingness
    :param info: Imputation Quality Confidence Score
    :return:
    '''
    out_dir = 'ukb_{}_maf{}_hwe{}_mind{}_geno{}_info{}'.format(imp_or_cal, maf, hwe, mind, geno, info)
    out_dir = out_dir.replace('.', '_').replace('e-', 'e')
    return out_dir


def create_subject_filter_list(genetic_sex_df, sex_df, aneuploidy_df, ethnicity_df, kinship_df):
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









