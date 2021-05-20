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
        parser.add_argument('--phesant-covar-csv-list', dest='phesant_covar_csv_list', action='append', nargs='+', default=None)
        parser.add_argument('--phesant-data-csv-list', dest='phesant_data_csv_list', default=None)
        parser.add_argument('--phesant-variablelistfile', dest='phesant_variablelistfile', action='append', nargs='+', default=None)
        parser.add_argument('--phesant-visits', dest='phesant_visits', default='0;1;2;3')
        parser.add_argument('--phesant-phenofile', dest='phesant_phenofile', action='append', nargs='+', default=None)
        parser.add_argument('--quiet', '-q', action='store_true', default=False)
        parser.add_argument('--noisy', '-n', action='count', default=0)

        args = parser.parse_args()
        return args


def main(argv=None):

    args = create_parser()
    # args = create_parser('yale')

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    configLogging(args)

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
    get_ukbutils(util_dir=os.path.join(config_json['base_dir'], 'data/ukb/raw'))

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
        prep_data_for_phesant(config_json, args)

    if 'run_phesant' in args.stage:
        run_phesant(config_json, args)

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
    csv_files     = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit0_regress.csv',
                     '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit2_regress.csv']

    table_types = dict(zip(['str','int','real','datetime'], ['VARCHAR', 'INTEGER', 'REAL', 'REAL']))

    ukb_meta = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv'))
    ukb_csv = csv_files[0]
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

        ukb_df.columns




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


def prep_data_for_phesant(config_json, args):
    '''Concat all bulk data into a single csv'''
    log.info('RUNNING: STAGE = PREP_DATA_FOR_PHESANT')

    # ukb showcase metadata important for aligning to age variables
    ukb_field_df = pd.read_table(os.path.join(config_json['repo_dir'], 'ref_files/field.txt'))

    # for now, just examine variables that were collected in person
    # age values need to be calculated by hand for online collections
    # see here: https://biobank.ndph.ox.ac.uk/showcase/schema.cgi?id=9
    inperson_fields = ukb_field_df['field_id'][ukb_field_df['instance_id'] == 2]
    grep_inperson_fields = ['x{}_'.format(x) for x in inperson_fields]
    # add imaging IDPs
    grep_inperson_fields = ['x25752_', 'x25750_', 'x25751_', 'x25753_', 'x25754_', 'x25755_']

    # use these ukb fields as covariates for phesant regressions
    phesant_covariates = ['21003',
                          '31',
                          '54',
                          '22000',
                          '22009',
                          '21003']

    # format them with 'x' prepended
    phesant_cols = ['xeid'] + \
                   ['x{}_'.format(covar) for covar in phesant_covariates]

    # read previously generated covariate csv files -- loop in case phesant covariates are split across multiple buckets
    log.debug('Found {} covariate files to read'.format(len(args.phesant_covar_csv_list)))
    log.debug('Found {} covariate files to read'.format(args.phesant_covar_csv_list))

    covar_list = []
    for covar_path in args.phesant_covar_csv_list:
        log.debug('Reading csv: {}'.format(covar_path[0]))
        covar_list.append(pd.read_csv(covar_path[0], low_memory=False))

    # concat all covar DFs and rename columns
    ukb_covar_df = pd.concat(covar_list)
    ukb_covar_df.columns = ['x{}'.format(x) for x in ukb_covar_df.columns]
    ukb_covar_df.columns = ukb_covar_df.columns.str.replace('-', '_').str.replace('[.]', '_')
    ukb_covar_cols = [x for x in list(ukb_covar_df.columns) if x != 'xeid']

    # read the big ukb dataframes
    log.debug('Found {} data files to read'.format(len(args.phesant_data_csv_list)))
    #for covar_path in args.phesant_covar_csv_list:
    ukb_csv = args.phesant_data_csv_list

    # read csv header
    hdr = pd.read_csv(ukb_csv, nrows=1)
    hdr.columns = ['x{}'.format(x) for x in hdr.columns.str.replace('-', '_').str.replace('[.]', '_')]
    hdr = hdr[hdr.columns[np.where(hdr.loc[0] != '#ukbgene')]]
    visitcol_dict = {}
    for visit in ['0','1','2','3']:
        visitcol_dict[visit] = {}
        this_visit_cols = hdr.columns[hdr.columns.str.contains('_{}_'.format(visit))]
        # TODO: add routines for online collections
        # only keep fields that were collected in-person
        notinperson_cols = [x for x in this_visit_cols if '{}_'.format(x.split('_')[0]) not in grep_inperson_fields]
        inperson_cols    = [x for x in this_visit_cols if '{}_'.format(x.split('_')[0]) in grep_inperson_fields]
        visitcol_dict[visit]['regress'] = ['xeid'] + inperson_cols
        visitcol_dict[visit]['process'] = ['xeid'] + notinperson_cols

    #ukb_csv = args.phesant_data_csv_list[0]
    step = 10000
    reader = pd.read_csv(ukb_csv, chunksize=step, low_memory=False, encoding="ISO-8859-1",)
    for i, ukb_df in enumerate(reader):
        log.debug('Reading Chunk-size: {}, Chunk: {}'.format(step, i))
        ukb_df.columns = ['x{}'.format(x) for x in ukb_df.columns]
        ukb_df.columns = ukb_df.columns.str.replace('-', '_').str.replace('[.]', '_')

        for visit in ['0','1','2','3']:
            column_info = visitcol_dict[visit]
            for type in ['regress', 'process']:
                # name of file to be written for phesant
                columns_to_write = column_info[type]
                phesant_write = ukb_csv.replace('.csv', '_phesant_visit{}_{}.csv'.format(visit, type))
                if i == 0: log.debug(phesant_write)

                # all columns for this visit (format ID_VISIT_INSTANCE
                vis_cols = list(ukb_df.columns[ukb_df.columns.str.contains('_{}_'.format(visit))])
                vis_df   = ukb_df[['xeid'] + vis_cols]

                if vis_df.shape[1] == 1: continue

                # remove columns in the covariate dataframe to prevent doubling up
                noncovar_cols = ['xeid'] + list(vis_df.columns.difference(ukb_covar_df.columns))

                # merge covariates with data
                merge_df   = ukb_covar_df.merge(vis_df[noncovar_cols], on='xeid')
                #visit_cols = visitcol_dict[visit]

                # rows, not including id columns
                check_cols = [x for x in columns_to_write if x != 'xeid']

                # get rid of subjects with all nans (e.g. non-imaging subjects on visit 2)
                keep_rows     = merge_df[check_cols].isna().sum(1) != len(check_cols)
                final_col_arr = ['xeid'] + pd.Series(list(set(check_cols + ukb_covar_cols))).sort_values().to_list()
                write_df      = merge_df[final_col_arr]
                write_df      = write_df.loc[keep_rows]
                # drop problematic columns that confuse phesant with commend characters
                #merge_df = merge_df[merge_df.columns[np.where(merge_df.loc[0] != '#ukbgene')]]
                if i == 0:
                    write_df.to_csv(phesant_write, index=None)
                else:
                    write_df.to_csv(phesant_write, mode='a', header=None, index=None)



def run_phesant(config_json, args):
    log.info('RUNNING: STAGE = RUN_PHESANT')
    log.debug(args)

    write_dir = os.path.join(config_json['base_dir'], 'slurm')

   # age_var     = args.age_var #'x21003_2_0'
    visit_arr = args.phesant_visits.split(';')

    # set up filepaths
    phesant_dir = os.path.join(config_json['repo_dir'], 'external/PHESANT/WAS')
    resDir      = os.path.join(config_json['base_dir'], 'data/ukb/phesant/')
    variablelistfile = os.path.join(config_json['repo_dir'], 'ref_files/outcome-info.tsv')
    datacodingfile   = os.path.join(config_json['repo_dir'], 'ref_files/data-coding-ordinal-info.txt')

    for visit in visit_arr:
        age_var = 'x21003_{}_0'.format(visit)
        for pheno_grep in args.phesant_phenofile:
            phenofile  = pheno_grep[0].replace('*', visit)
            phenoname  = phenofile.split('/')[-1].replace('.csv','')
            phesant_cmd = f'''
cd {phesant_dir}
Rscript phenomeScan.r \\
        --phenofile={phenofile} \\
        --variablelistfile={variablelistfile} \\
        --datacodingfile={datacodingfile} \\
        --traitofinterest="{age_var}" \\
        --resDir={resDir} \\
        --userId="xeid" \\
        --sensitivity \\
        --ageistraitofinterest \\
        --genetic="FALSE" \\
        --visit={visit} \\
        --file_prepend={phenoname} \\
        --save'''

            # write executable to file
            if args.slurm == True:
                write_file = os.path.join(write_dir, 'phesant_{}'.format(phenoname))
                utilities.submit_slurm(
                    utilities.write_slurm(slurm_file=write_file,
                                          partition=args.slurm_partition, cmd=phesant_cmd,
                                          jobName=phenoname, stime='6:00:00', nthreads=12))
            else:
                write_file = os.path.join(write_dir, 'phesant_{}.bash'.format(phenoname))
                utilities.write_cmd(cmd=phesant_cmd, write_file=write_file)


def snp_preprocess(config_json, args):
    args.maf = 0.01
    args.hwe = 1e-6
    args.mind = 0.02
    args.geno = 0.02
    args.info = 0.60

    # point to directories
    genotyped_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/genotyped')
    imputed_dir   = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed')

    # create descriptive folder name
    imp_dirname  = genetics.snp_dirname(imp_or_cal='imp', maf=args.maf, hwe=args.hwe, mind=args.mind, geno=args.geno, info=args.info)
    snp_dirname  = genetics.snp_dirname(imp_or_cal='cal', maf=args.maf, hwe=args.hwe, mind=args.mind, geno=args.geno, info=args.info)

    # create output directory
    geno_proc_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/genotyped_processed')
    if not os.path.exists(geno_proc_dir):
        os.mkdir(geno_proc_dir)
    imp_proc_dir = os.path.join(config_json['base_dir'], 'data/ukb/genetic/imputed_processed')
    if not os.path.exists(imp_proc_dir):
        os.mkdir(imp_proc_dir)




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



def configLogging(args):
    """Configures ``funpack`` logging.

    :arg args: ``argparse.Namespace`` object containing parsed command line
               arguments.
    """
    class LogHandler(logging.StreamHandler):
        def emit(self, record):
            levelno = record.levelno
            if   levelno >= logging.WARNING:  colour = '\x1b[31;1m'
            elif levelno >= logging.INFO:     colour = '\x1b[39;1m'
            elif levelno >= logging.DEBUG:    colour = '\x1b[90;1m'
            else:                             colour = ''
            # Reset terminal attributes
            # after each message.
            record.msg = '{}{}\x1b[0m'.format(colour, record.msg)
            return super(LogHandler, self).emit(record)

    fmt    = logging.Formatter('%(asctime)s '
                               '%(levelname)8.8s '
                               '%(filename)20.20s '
                               '%(lineno)4d: '
                               '%(funcName)-15.15s - '
                               '%(message)s',
                               '%H:%M:%S')
    handler = LogHandler()
    handler.setFormatter(fmt)
    log.addHandler(handler)
    log.propagate = False



if __name__ == '__main__':
    sys.exit(main())






