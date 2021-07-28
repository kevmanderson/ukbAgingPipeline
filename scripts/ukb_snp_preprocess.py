#!/bin/python

import os
import json
import shutil
import subprocess
import argparse
import sqlite3
import glob
import logging
import datetime
import numpy as np
import pandas as pd
import sqlite3
import utilities.utilities as utilities

log = logging.getLogger('ukb')

def main(config_json, args):
    '''

    :param config_json:
    :param args:
    :return:
    '''

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    log.info('RUNNING: STAGE = SNP PREPROCESSING')

    parser = argparse.ArgumentParser()
    parser.add_argument('--ukb_csv', dest='ukb_csv', required=True, default=None)
    parser.add_argument('--bgens', dest='bgens', required=True, default=None)
    parser.add_argument('--samples', dest='samples', required=True, default=None)
    parser.add_argument('--info_files', dest='info_files', required=True, default=None)
    parser.add_argument('--maf', dest='maf', required=True, default=None)
    parser.add_argument('--hwe', dest='hwe', required=True, default=None)
    parser.add_argument('--mind', dest='mind', required=True, default=None)
    parser.add_argument('--geno', dest='geno', required=True, default=None)
    parser.add_argument('--info', dest='info', required=True, default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True, default=None)
    parser.add_argument('--slurm', dest='slurm', action='store_true', default=False)
    parser.add_argument('--slurm_dir', dest='slurm_dir', required=True, default=None)
    parser.add_argument('--plink2', dest='out_dir', required=True, default=None)

    parser = argparse.ArgumentParser()
    args   = parser.parse_args()

    tmp = False
    if tmp == True:
        args.ukb_csv = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/phesant_visit0.csv'
        args.bgens = '/gpfs/milgram/data/UKB/ukb_snp/ukb_imp_chr*_v3.bgen'
        args.samples = '/gpfs/milgram/data/UKB/ukb_snp/ukb25163_imp_chr*_v3_s487324.sample'
        args.info_files = '/gpfs/milgram/data/UKB/ukb_snp/ukb_mfi_chr*_v3.txt'
        args.maf = 0.01
        args.hwe = 1e-6
        args.mind = 0.1
        args.geno = 0.1
        args.info = 0.8
        args.slurm = True
        args.slurm_dir = '/gpfs/milgram/scratch/holmes/simons_aging/slurm'
        args.plink2  = '/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2'
        args.out_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/genetic/imputed_processed'

    # list of autosomal bgens
    bgen_files = [args.bgens.replace('*',str(i)) for i in range(1,23)]

    # create descriptive folder names
    imp_dirname = snp_dirname(imp_or_cal='imp',
                               maf=args.maf,
                               hwe=args.hwe,
                               mind=args.mind,
                               geno=args.geno,
                               info=args.info)

    # create output directory
    utilities.make_dir(os.path.join(args.out_dir, imp_dirname))

    # identify subjects for polygenic processing
    cohort_filter(ukb_csv=args.ukb_csv, args)

    # make plink-formatted subject list for those passing initial QC filters
    plink_subfile = cohort_filter(ukb_csv=args.ukb_csv,
                                  rm_sex_mismatch=True,
                                  rm_sex_aneuploidy=True,
                                  rm_non_eur=True,
                                  rm_related=True,
                                  rm_het_outliers=True,
                                  out_dir=args.out_dir)


    sub_keep_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/plink_genetic_subjects.txt')
    plink2 = '/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2'

    for chrom in range(1,23):
        log.info('Working on CHR={}'.format(chrom))
        #fname      = imp_dirname.replace('_imp_', '_imp_chr{}_'.format(chrom))
        fname      = '{}_chr{}'.format(imp_dirname, chrom)
        bgen_out   = os.path.join(args.out_dir, imp_dirname, fname)

        # chromosomal bgen/sample
        cur_bgen   = args.bgens.replace('*', str(chrom))
        cur_sample = args.samples.replace('*', str(chrom))
        cur_info   = args.info_files.replace('*', str(chrom))

        log.info('Reading MFI (INFO score) for chrom {}'.format(chrom))
        info_df = pd.read_table(cur_info, sep='\t', header=None)
        info_df.columns = ['id', 'rsid', 'bp', 'a1', 'a2', 'maf', 'f', 'info']
        keep_info_df = info_df.loc[info_df['info'] > args.info]

        # write plink range file for SNPs to retain
        info_write_df = keep_info_df[['bp', 'bp']]
        info_write_df.insert(0, 'chr', chrom)
        snp_extract_file = os.path.join(args.out_dir, imp_dirname, 'plink_snp_range_chr{}.txt'.format(chrom))
        info_write_df.to_csv(snp_extract_file, index=None, sep='\t', header=None)

        # create plink command for chromsomal QC
        plink_cmd = f'''{args.plink2} \\
         --bgen {cur_bgen} ref-first \\
         --sample {cur_sample} \\
         --max-alleles 2 \\
         --maf {args.maf} \\
         --hwe {args.hwe} \\
         --mind {args.mind} \\
         --geno {args.geno} \\
         -export bgen-1.2 'bits=8' \\
         --keep {plink_subfile} \\
         --extract 'range' {snp_extract_file} \\
         --out {bgen_out}
         '''

        # submit this to slurm b/c parallelization is cool
        if args.slurm == True:
            slurm_file = os.path.join(args.slurm_dir, 'snp_preproc_chrom{}'.format(chrom))
            slurm_path = utilities.write_slurm(slurm_file=slurm_file,
                                               partition='psych_day',
                                               cmd=plink_cmd,
                                               jobName='plink_chr{}'.format(chrom),
                                               stime='6:00:00',
                                               n_gpu=None,
                                               nthreads=20,
                                               mem='140GB')
            log.info(slurm_path)
            job_id = utilities.submit_slurm(slurm_path)
        # or just write the plink command to a text file
        else:
            utilities.write_cmd(cmd=plink_cmd, write_file=slurm_file)


def cohort_filter(ukb_csv,
                  rm_sex_mismatch,
                  rm_sex_aneuploidy,
                  rm_non_eur,
                  rm_related,
                  rm_het_outliers,
                  out_dir):
    '''
    :cvar
    '''
    log.info('Identifying UKB subjects for genetic analyses')

    # fields required for cohort filtering
    # -----------
    genotype_qc_fields = ['31', '22019', '22001', '22021', '22006', '22004', '22005', '22020', '22009', '22027']
    geno_qc_substrings = ['x{}_'.format(x) for x in genotype_qc_fields]

    # path to previously generated dataframe (combines all UKB buckets)
    log.info('Reading genetic covariates from: {}'.format(ukb_csv))

    # read header and make list of variables to pull
    with open(ukb_csv, 'r') as f:
        hdr = pd.Series(f.readline().replace('\n', '').split(','))

    # indices/column-names for the genetic data to extract
    # -----------
    idx_lists    = [list(np.where(hdr.str.contains(x))[0]) for x in geno_qc_substrings]
    flat_list    = [item for sublist in idx_lists for item in sublist]
    pull_columns = hdr[flat_list].tolist()

    # read genetic quality control information from unified visit0 csv file
    # speeds up the read to only retain a subset of the data
    # -----------
    log.info('Start read...')
    gene_qc_dfs = []
    reader = pd.read_csv(ukb_csv, chunksize=10000, low_memory=False, encoding="ISO-8859-1",)
    for i, ukb_df in enumerate(reader):
        print(i)
        gene_qc_dfs.append(ukb_df[['xeid'] + pull_columns])
    gene_df = pd.concat(gene_qc_dfs)
    log.info('Done.')


    # running array of subjects to remove
    rm_subs = []


    # GENETIC / REPORTED SEX mismatch
    # -----------
    if rm_sex_mismatch:
        rm_sexmismatch = gene_df.loc[gene_df['x22001_0_0'] != gene_df['x31_0_0'], 'xeid'].tolist()
        log.info('Reported-vs-Genetic sex mismatch: n={}'.format(len(rm_sexmismatch)))
        rm_subs = rm_subs + rm_sexmismatch

    # SEX ANEUPLOIDY
    # -----------
    if rm_sex_aneuploidy:
        rm_aneuploidy = gene_df.loc[np.where(gene_df['x22019_0_0'] == 1), 'xeid'].tolist()
        log.info('Sex aneuploidy: n={}'.format(len(rm_aneuploidy)))
        rm_subs = rm_subs + rm_aneuploidy

    # EUROPEAN
    # -----------
    if rm_non_eur:
        rm_ethnicity = gene_df.loc[gene_df['x22006_0_0'] != 1, 'xeid'].tolist()
        log.info('Not EUR ancestry: n={}'.format(len(rm_ethnicity)))
        rm_subs = rm_subs + rm_ethnicity

    # KINSHIP
    # -----------
    #kin_df = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_rel_a25163_s488225.dat'), delim_whitespace=True)
    #kin_df['idx'] = kin_df.index
    # third degree relatives
    #cutoff = 0.0884
    #kin_df = kin_df.loc[kin_df.Kinship > cutoff]
    #connections = kin_df[['ID1','ID2']].melt()['value'].value_counts()

    # for now, use the UKB provided "used.in.pc.calculation', which is a maximal set of unrelated subjects
    if rm_related:
        rm_related = gene_df.loc[np.where(gene_df['x22020_0_0'] != 1)[0], 'xeid'].tolist()
        log.info('Not used in genetic PC calc (relatedness): n={}'.format(len(rm_related)))
        rm_subs = rm_subs + rm_related

    # HETEROZYGOSITY / MISSINGNESS OUTLIERS
    # -----------
    if rm_het_outliers:
        rm_outliers = gene_df.loc[np.where(gene_df['x22027_0_0'] == 1)[0], 'xeid'].tolist()
        log.info('Het/Missingness outliers: n={}'.format(len(rm_outliers)))
        rm_subs = rm_subs + rm_outliers

    rm_subs = list(set(rm_subs))
    log.info('UNIQUE subs to remove: n={}'.format(len(rm_subs)))

    # subjects for genetic analyses
    # -----------
    subs_for_genetics = gene_df.loc[~gene_df.xeid.isin(rm_subs), 'xeid'].tolist()
    plink_df = pd.DataFrame({'IID': subs_for_genetics, 'FID': subs_for_genetics})
    log.info('UKB subjects remaining: n={}'.format(plink_df.shape[0]))

    write_file = os.path.join(out_dir, 'genetic_subjects.txt')
    pd.Series(subs_for_genetics).to_csv(write_file, index=None, header=None)

    write_file = os.path.join(out_dir, 'plink_genetic_subjects.txt')
    log.info('PLINK formatted subject keep file: {}'.format(write_file))
    plink_df.to_csv(write_file, index=None, sep='\t')

    return write_file



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

def info_filter(config_json, args):
    '''
    :cvar
    '''

    log.info('Filtering imputed variants with INFO <= {}'.format(args.info))
    imputed_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed')
    qctool = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/qctool/qctool'

    for chrom in range(1,23):
        info_file = os.path.join(imputed_dir, 'ukb_mfi_chr{}_v3.txt'.format(chrom))
        log.info('Reading: {}'.format(info_file))

        info_df   = pd.read_table(info_file, header=None)
        info_df.columns = ['id', 'rsid', 'bp', 'ref', 'alt', 'maf', 'fill', 'info']
        log.info('Total variants: {}'.format(info_df.shape[0]))

        info_keep_df  = info_df.loc[info_df['info'] > args.info]
        keep_var_df   = pd.DataFrame({'chr1': chrom, 'chr2': chrom, 'bp': info_keep_df.bp.astype(str)})
        log.info('Passing variants: {}'.format(info_keep_df.shape[0]))

        info_file = os.path.join(imputed_dir, 'ukb_mfi_chr{}_v3_plink_keep.txt'.format(chrom))
        log.info('Writing: {}'.format(info_file))
        keep_var_df.to_csv(info_file, index=None, header=None, sep='\t')




if __name__ == "__main__":
    main()






