#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# this is the primary entrypoint for the 'buckner-lab-ukb-pipeline' repo
#
"""This module contains the ```buckner-lab-ukb-pipeline``` entrypoint"""

import os
import sys
import json
import warnings
import argparse
import numpy as np
import logging
import datetime
import calendar
import pandas as pd
import glob
import wget
import tarfile
from functools import reduce

import phesant
import download
import genetics
import utilities.utilities as utilities
from utilities.utilities import write_slurm, submit_slurm, make_dir, make_symlink
from download import decrypt_ukb_data, get_ukbutils, read_restbulk_and_write
log = logging.getLogger('ukb')

def create_parser(interactive=None):
    example_text = '''example: 
     python test.py -t template/test.py'''

    # tmp manual definition of parser arguments
    if interactive == 'yale':
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.base_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging'
        args.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json'
        args.bulk_field = [['rfmri_full_25:25750'], ['rfmri_full_100:25751'], ['rfmri_part_25:25752'], ['rfmri_part_100:25753'], ['rfmri_rsfa_25:25754'], ['rfmri_rsfa_100:25755']]
        args.make_bulk_list = True
        args.download_bulk_data = False
        args.slurm = True
        args.stage = 'convert'
        args.slurm_partition = 'short'
        args.convert_all_fields = False
        args.log_file = None
        args.quiet = False
        args.noisy = 0
        args.phesant_data_csv_list = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_2.csv'
        args.phesant_covar_csv_list = [['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_covars.csv']]
        args.phesant_data_csv_list = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.csv',
                                      '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.csv']
        args.singularity_container = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/simons_ukb_aging_pipeline'
        args.phesant_phenofile = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit*_regress.csv',
                                  '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit*_process.csv',
                                  '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb43410_phesant_visit*_regress.csv',
                                  '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb43410_phesant_visit*_process.csv']
        args.phesant_csv = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_3.csv']
        args.phesant_visits = '0;1;2;3'
        return args

    # tmp manual definition of parser arguments
    elif interactive == 'harvard':
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.base_dir = '/ncf/sba01/simons_aging'
        args.config = '/ncf/sba01/ukbAgingPipeline/harvard_config.json'
        args.bulk_field = [['rfmri_full_25:25750'], ['rfmri_full_100:25751'], ['rfmri_part_25:25752'], ['rfmri_part_100:25753'], ['rfmri_rsfa_25:25754'], ['rfmri_rsfa_100:25755']]
        args.make_bulk_list = True
        args.download_bulk_data = False
        args.slurm = True
        args.stage = 'convert'
        args.slurm_partition = 'short'
        args.convert_all_fields = False
        args.log_file = None
        args.quiet = False
        args.noisy = 0
        args.phesant_data_csv_list = '/ncf/sba01/simons_aging/data/ukb/raw/bulk_25750_2.csv'
        args.phesant_visits = '0;1;2;3'
        return args

    else:
        parser = argparse.ArgumentParser()
        parser.add_argument('--config', '-c', dest='config', required=True)
        parser.add_argument('--stage', dest='stage', required=False)
        parser.add_argument('--bulk-field', '-b', dest='bulk_field', action='append', nargs='+', required=False)
        parser.add_argument('--make-bulk-list', dest='make_bulk_list', action='store_true', default=False)
        parser.add_argument('--download-bulk-data', dest='download_bulk_data', action='store_true', default=False)
        parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
        parser.add_argument('--slurm_partition', '-p', dest='slurm_partition', required=False, default='short')
        parser.add_argument('--singularity_container', dest='singularity_container', required=False, default=False)
        parser.add_argument('--stages', dest='stages', nargs='+', required=False, default=False)
        parser.add_argument('--convert_all_fields', action='store_true', required=False, default=False)
        parser.add_argument('--log_file', dest='log_file', default=None)
        parser.add_argument('--maf', dest='maf', type=float, default=None)
        parser.add_argument('--info', dest='info', type=float, default=None)
        parser.add_argument('--mind', dest='mind', type=float, default=None)
        parser.add_argument('--geno', dest='geno', type=float, default=None)
        parser.add_argument('--hwe', dest='hwe', type=float, default=None)
        parser.add_argument('--phesant-covar-csv-list', dest='phesant_covar_csv_list', action='append', nargs='+', default=None)
        parser.add_argument('--phesant-data-csv-list', dest='phesant_data_csv_list', default=None)
        parser.add_argument('--phesant-csv', dest='phesant_csv', default=None)
        parser.add_argument('--phesant-variablelistfile', dest='phesant_variablelistfile', action='append', nargs='+', default=None)
        parser.add_argument('--phesant-visits', dest='phesant_visits', default='0;1;2;3')
        parser.add_argument('--phesant-phenofile', dest='phesant_phenofile', action='append', nargs='+', default=None)
        parser.add_argument('--quiet', '-q', action='store_true', default=False)
        parser.add_argument('--noisy', '-n', action='count', default=0)
        parser.add_argument('--use-enc-ukb', dest='use_enc_ukb', default=None)
        args = parser.parse_args()
        return args


def main(argv=None):

    args = create_parser()
    # args = create_parser('yale')
    # args = create_parser('harvard')

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    utilities.configLogging(args)

    log.info('buckner-lab-ukb-pipeline')
    log.info('Date: %s (%s)', date.today(), calendar.day_name[date.weekday()])

    # read configuration file
    log.debug('Reading config file here: {}'.format(args.config))
    with open(args.config, 'r') as f:
        config_json = json.load(f)[0]
    slurm_dir = os.path.join(config_json['base_dir'], 'slurm')

    # create project directory structure if necessary
    utilities.create_directories(root_dir=config_json['base_dir'])

    # download UKB utilities if needed
    get_ukbutils(util_dir=os.path.join(config_json['base_dir'], 'data/ukb/external'))

    # decrypt the encoded ukb data
    if 'decrypt' == args.stage:
        download.decrypt_stage(config_json, args)

    # convert encoded data files
    if 'convert' in args.stage:
        download.convert_stage(config_json, args)

    # make list of subjects to download for bulk
    if 'make_bulk_list' in args.stage:
        download.make_bulk_subj_list(config_json, args)

    # download bulk field
    if 'download_bulk' in args.stage:
        download.download_bulk_stage(config_json, args)

    # compile already downloaded bulk data to csv
    if 'bulk_to_csv' in args.stage:
        download.bulk_to_csv(config_json, args)

    if 'prep_data_for_phesant' in args.stage:
        phesant.prep_data_for_phesant(config_json, args)

    if 'run_phesant' in args.stage:
        phesant.run_phesant(config_json, args)

    if 'ukb_sql' in args.stage:
        ukb_sql(config_json, args)

    if 'download_genetic' in args.stage:
        download_genetic(config_json, args)

    if 'snp_preprocess' in args.stage:
        snp_preprocess(config_json, args)

    if 'make_metadata' in args.stage:
        make_metadata(config_json, args)



def field_to_tables(table_types, ukb_meta, ukb_cols):

    ukb_cols_short = ['{}_'.format(x.split('_')[0]) for x in ukb_cols]
    ukb_col_map = dict(zip(ukb_cols_short, ukb_cols))
    ukb_cols_short = [x for x in ukb_cols_short if x != 'xeid_']

    table_field_map = {}
    table_map = dict(zip(['str', 'int', 'real', 'datetime'], [['Categorical multiple','Categorical single'],
                                                              ['Integer'],
                                                              ['Continuous'],
                                                              ['Date']]))
    for table_name, sql_type in table_types.items():
        print('{}/{}'.format(table_name, sql_type))
        table_fields = ukb_meta['fieldid'][ukb_meta['valuetype'].isin(table_map[table_name])].tolist()
        table_fields_fmt = ['x{}_'.format(x) for x in table_fields]
        table_fields_fmt = list(set(table_fields_fmt).intersection(set(ukb_cols_short)))
        save_fields = [ukb_col_map[x] for x in table_fields_fmt]
        table_field_map[table_name] = save_fields

    return table_field_map


def insert_data_to_sql(ukb_df, table_name, table_fields):
    table_cols = ['x{}_'.format(x) for x in table_fields]
    curr_tab_fields = set(ukb_df.columns.to_list()).intersection(set(table_fields))
    trips = chunk[['eid'] + list(curr_tab_fields)].melt(id_vars='eid', value_vars=curr_tab_fields)
    trips = trips[trips['value'].notnull()]


def ukb_sql(config_json, args):

    sql_dir       = os.path.join(config_json['base_dir'], 'data/ukb/sql')
    repo_dir      = config_json['repo_dir']

    ukb_meta      = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv'))
    ordinal_codes = pd.read_csv(os.path.join(repo_dir, 'ref_files/data-coding-ordinal-info.txt'))
    codings       = pd.read_csv(os.path.join(repo_dir, 'ref_files/codings.csv'))

    col_list = []
    csv_files = glob.glob(os.path.join(config_json['base_dir'], 'data/ukb/phesant', '*-data-*30.txt'))
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        col_list = col_list + df.columns.tolist()
    flat_list = [item for sublist in col_list for item in sublist]

    table_types = dict(zip(['str','int','real','datetime'], ['VARCHAR', 'INTEGER', 'REAL', 'REAL']))

    ukb_meta = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv'))
    ukb_csv  = csv_files[0]
    sql_path = os.path.join(sql_dir, 'ukb_db.sqlite')
    step = 10000

    hdr = pd.read_csv(ukb_csv, nrows=1, encoding="ISO-8859-1")
    table_cols = field_to_tables(table_types=table_types, ukb_meta=ukb_meta, ukb_cols=hdr.columns.tolist())

    reader = pd.read_csv(ukb_csv, chunksize=step, low_memory=False, encoding="ISO-8859-1",)
    for i, ukb_df in enumerate(reader):
        log.debug('Reading Chunk-size: {}, Chunk: {}'.format(step, i))
        chunk_cols = set(ukb_df.columns) - set('xeid')
        for table_name, sql_type in table_types.items():
            table_fields = table_cols[table_name]
            x = insert_data_to_sql(ukb_df, table_name, table_fields)
        field_desc.to_sql("field_desc", con, if_exists='replace', index=False)

        #ukb_df.columns


def make_metadata(config_json, args):

    repo_dir = config_json['repo_dir']

    showcase      = pd.read_csv(os.path.join(repo_dir, 'ref_files/showcase.csv'))
    codings       = pd.read_csv(os.path.join(repo_dir, 'ref_files/codings.csv'))
    outcome_info  = pd.read_table(os.path.join(repo_dir, 'ref_files/outcome-info.tsv'))
    ordinal_codes = pd.read_csv(os.path.join(repo_dir, 'ref_files/data-coding-ordinal-info.txt'))
    regr_filter   = pd.read_csv(os.path.join(repo_dir, 'ref_files/ukb_allow_regression.csv'))

    fields        = pd.read_table(os.path.join(repo_dir, 'ref_files/field.txt'))
    fields.columns = fields.columns.str.replace('field_id','fieldid')

    # outcome_info comes from PHESANT. use this as the base reference dataframe, and add from there
    outcome_info.columns = outcome_info.columns.str.lower()
    outcome_cols = [
        'fieldid',
        'excluded',
        'cat_mult_indicator_fields',
        'cat_single_to_cat_mult',
        'date_convert',
        'valuetype'
    ]

    # showcase data from ukb website; https://biobank.ndph.ox.ac.uk/showcase/schema.cgi?id=1
    showcase.columns = showcase.columns.str.lower()
    showcase_cols = ['fieldid',
                     'path',
                     'category',
                     'field',
                     'participants',
                     'items',
                     'units',
                     'sexed',
                     'coding',
                     'notes',
                     'link',
                     ]

    showcase_merge = showcase[showcase_cols]
    outcome_info_merge = outcome_info[outcome_cols]

    # merge all the data
    meta_df = fields.merge(showcase_merge, on='fieldid')
    meta_df = meta_df.merge(outcome_info_merge, on='fieldid')
    meta_df = meta_df.merge(regr_filter[['fieldid','allow_regression']], on='fieldid')

    write_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv')
    meta_df.to_csv(write_file, index=None)
    return meta_df




def snp_preprocess(config_json, args):
    '''

    :param config_json:
    :param args:
    :return:
    '''
    #args.maf = 0.01
    #args.hwe = 1e-6
    #args.mind = 0.02
    #args.geno = 0.02
    #args.info = 0.60

    # point to genotype data directories
    # --------------------
    genotyped_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/genotyped')
    imputed_dir   = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed')
    log.info('GENOTYPED dir : {}'.format(genotyped_dir))
    log.info('IMPUTED dir : {}'.format(genotyped_dir))

    # create descriptive folder names
    imp_dirname  = genetics.snp_dirname(imp_or_cal='imp', maf=args.maf, hwe=args.hwe, mind=args.mind, geno=args.geno, info=args.info)
    snp_dirname  = genetics.snp_dirname(imp_or_cal='cal', maf=args.maf, hwe=args.hwe, mind=args.mind, geno=args.geno, info=args.info)

    # create output directories
    # --------------------
    geno_proc_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/genotyped_processed')
    utilities.make_dir(geno_proc_dir)
    utilities.make_dir(os.path.join(geno_proc_dir, snp_dirname))

    imp_proc_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed_processed')
    utilities.make_dir(imp_proc_dir)
    utilities.make_dir(os.path.join(imp_proc_dir, imp_dirname))

    # identify high quality imputed SNPs
    #info_filter(config_json, args)
    #cohort_filter(config_json, args)


    sub_keep_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/plink_genetic_subjectsar.txt')
    plink2 = '/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2'
    for chrom in range(1,23):
        fname      = imp_dirname.replace('_imp_', '_imp_chr{}_'.format(chrom))
        bgen_out   = os.path.join(imp_proc_dir, imp_dirname, fname)

        cur_bgen   = '/gpfs/milgram/data/UKB/ukb_snp/ukb_imp_chr{}_v3.bgen'.format(chrom)
        cur_sample = '/gpfs/milgram/data/UKB/ukb_snp/ukb25163_imp_chr{}_v3_s487324.sample'.format(chrom)
        snp_keep_file  = os.path.join(imputed_dir, 'ukb_mfi_chr{}_v3_plink_keep.txt'.format(chrom))
        sub_keep_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/plink_genetic_subjects.txt')

        plink_cmd = f'''{plink2} \\
         --bgen {cur_bgen} ref-first \\
         --sample {cur_sample} \\
         --max-alleles 2 \\
         --maf {args.maf} \\
         --hwe {args.hwe} \\
         --mind {args.mind} \\
         --geno {args.geno} \\
         -export bgen-1.2 'bits=8' \\
         --keep {sub_keep_file} \\
         --extract 'range' {snp_keep_file} \\
         --out {bgen_out}
         '''

        slurm_dir  = os.path.join(config_json['base_dir'], 'slurm')
        slurm_file = os.path.join(slurm_dir, 'snp_preproc_chrom{}'.format(chrom))
        if args.slurm == True:
            slurm_path = utilities.write_slurm(slurm_file, 'short', plink_cmd, 'plink_chr{}'.format(chrom), stime='6:00:00', n_gpu=None, nthreads=24, mem='180G')
            log.info(slurm_path)
            job_id = utilities.submit_slurm(slurm_path)

        else:
            utilities.write_cmd(cmd=plink_cmd, write_file=slurm_file)


    exit()


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


def cohort_filter(config_json, args):
    '''
    :cvar
    '''

    log.info('Identifying UKB subjects for genetic analyses')

    # fields required for genotype filter
    genotype_qc_fields = ['31', '22019', '22001', '22021', '22006', '22004', '22005', '22020', '22009', '22027']
    geno_qc_substrings = ['x{}_'.format(x) for x in genotype_qc_fields]

    # path to previously generated dataframe (combines all UKB buckets)
    vis0_path = os.path.join(config_json['base_dir'], 'data/ukb/raw/phesant_visit0.csv')
    log.info('Reading genetic covariates from: {}'.format(vis0_path))

    # read header and make list of variables to pull
    with open(vis0_path, 'r') as f:
        hdr = pd.Series(f.readline().replace('\n', '').split(','))

    idx_lists = [list(np.where(hdr.str.contains(x))[0]) for x in geno_qc_substrings]
    flat_list = [item for sublist in idx_lists for item in sublist]
    pull_columns = hdr[flat_list].tolist()

    # read genetic quality control information from unified visit0 dataframe
    log.info('Start read...')
    gene_qc_dfs = []
    reader = pd.read_csv(vis0_path, chunksize=10000, low_memory=False, encoding="ISO-8859-1",)
    for i, ukb_df in enumerate(reader):
        gene_qc_dfs.append(ukb_df[['xeid'] + pull_columns])
    gene_df = pd.concat(gene_qc_dfs)
    log.info('Done.')

    # GENETIC vs REPORTED SEX
    rm_sexmismatch = gene_df.loc[gene_df['x22001_0_0'] != gene_df['x31_0_0'], 'xeid'].tolist()
    log.info('Reported-vs-Genetic sex mismatch: n={}'.format(len(rm_sexmismatch)))


    # SEX ANEUPLOIDY
    rm_aneuploidy = gene_df.loc[np.where(gene_df['x22019_0_0'] == 1), 'xeid'].tolist()
    log.info('Sex aneuploidy: n={}'.format(len(rm_aneuploidy)))

    # EUROPEAN
    rm_ethnicity = gene_df.loc[gene_df['x22006_0_0'] != 1, 'xeid'].tolist()
    log.info('Not EUR ancestry: n={}'.format(len(rm_ethnicity)))

    # KINSHIP
    #kin_df = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_rel_a25163_s488225.dat'), delim_whitespace=True)
    #kin_df['idx'] = kin_df.index
    # third degree relatives
    #cutoff = 0.0884
    #kin_df = kin_df.loc[kin_df.Kinship > cutoff]
    #connections = kin_df[['ID1','ID2']].melt()['value'].value_counts()

    # for now, use the UKB provided "used.in.pc.calculation', which is a maximal set of unrelated subjects
    rm_related = gene_df.loc[np.where(gene_df['x22020_0_0'] != 1)[0], 'xeid'].tolist()
    log.info('Not used in genetic PC calc (relatedness): n={}'.format(len(rm_related)))

    # HETEROZYGOSITY / MISSINGNESS OUTLIERS
    rm_outliers = gene_df.loc[np.where(gene_df['x22027_0_0'] == 1)[0], 'xeid'].tolist()
    log.info('Het/Missingness outliers: n={}'.format(len(rm_outliers)))

    rm_subs = rm_sexmismatch + rm_aneuploidy + rm_ethnicity + rm_related + rm_outliers
    rm_subs = list(set(rm_subs))
    log.info('UNIQUE subs to remove: n={}'.format(len(rm_subs)))

    # subjects for genetic analyses
    subs_for_genetics = gene_df.loc[~gene_df.xeid.isin(rm_subs), 'xeid'].tolist()
    plink_df = pd.DataFrame({'IID': subs_for_genetics, 'FID': subs_for_genetics})
    log.info('UKB subjects remaining: n={}'.format(plink_df.shape[0]))

    write_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/genetic_subjects.txt')
    pd.Series(subs_for_genetics).to_csv(write_file, index=None, header=None)

    write_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/plink_genetic_subjects.txt')
    log.info('PLINK formatted subject keep file: {}'.format(write_file))
    plink_df.to_csv(write_file, index=None, sep='\t')



def download_genetic(config_json, args):

    gfetch = os.path.join(config_json['base_dir'], 'data/ukb/external/gfetch')

    # list of chromosomes to download
    chrom_list = [str(chr) for chr in range(1,27)]
    chrom_list = chrom_list + ['X', 'Y', 'XY', 'MT']

    # point to directories
    genotyped_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/genotyped')
    imputed_dir   = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed')

    # if the user specified already downloaded data
    if 'genotyped_data' in config_json.keys():
        geno_files = glob.glob('{}*'.format(config_json['genotyped_data']))
        for src_file in geno_files:
            print(geno_file)
            fname = src_file.split('/')[-1]
            dest_file = os.path.join(genotyped_dir, fname)
            make_symlink(src_file, dest_file)

    # if the user specified already downloaded data
    if 'imputed_data' in config_json.keys():
        bgen_files   = glob.glob('{}.bgen'.format(config_json['imputed_data'][0]['bgen']))
        sample_files = glob.glob('{}.sample'.format(config_json['imputed_data'][0]['sample']))
        for src_file in bgen_files + sample_files:
            fname = src_file.split('/')[-1]
            dest_file = os.path.join(imputed_dir, fname)
            make_symlink(src_file, dest_file)

    # copy the ukbgene utility
    orig_ukbgene = os.path.join(config_json['repo_dir'], 'external/ukbgene')
    shutil.copy(orig_ukbgene, os.path.join(genotyped_dir, 'ukbgene'))
    shutil.copy(orig_ukbgene, os.path.join(imputed_dir, 'ukbgene'))

    # SNP QC data
    wget.download('http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_snp_qc.txt', out=genotyped_dir)
    os.chmod(os.path.join(genotyped_dir, 'ukb_snp_qc.txt'), 755)

    # IMPUTED info
    wget.download('http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_imp_mfi.tgz', out=imputed_dir)
    os.chdir(imputed_dir)
    tar = tarfile.open(os.path.join(imputed_dir, 'ukb_imp_mfi.tgz'))
    tar.extractall()
    tar.close()

    # relatedness data
    #rel_fetch = f'''cd {genotyped_dir}\n'''
    #rel_fetch = f'''{rel_fetch}{gfetch} rel'''

    #./gfetch rel -aukb40501.key



def convert_cmd(ukbconv, raw_dir, enc_ukb, output_format, field_file=None, out=None):
    '''Create command to convert ukb data fields'''
    print(out)
    if out is not None:
        o_name = out
    else:
        o_name = enc_ukb.replace('.enc_ukb','')
    decrypt_cmd = f"cd {raw_dir}\n"
    if field_file is not None:
        decrypt_cmd = f"{decrypt_cmd}{ukbconv} {enc_ukb} {output_format} -o./{o_name} -i./{field_file}"
    else:
        decrypt_cmd = f"{decrypt_cmd}{ukbconv} {enc_ukb} {output_format} -o./{o_name}"
    return decrypt_cmd


if __name__ == '__main__':
    sys.exit(main())






