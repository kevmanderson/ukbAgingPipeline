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
import copy
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
import modules.ukb_preprocess as ukb_preprocess
import modules.ukb_genetics as ukb_genetics
import modules.ukb_phesant as ukb_phesant
import modules.ukb_sql as ukb_sql

#import download
#import genetics
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
        args.maf = 0.01
        args.info = 0.80
        args.mind = 0.10
        args.geno = 0.10
        args.hwe = 1e6
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
        ukb_preprocess.decrypt_stage(config_json, args)

    # convert encoded data files
    if 'convert' in args.stage:
        ukb_preprocess.convert_stage(config_json, args)

    # make list of subjects to download for bulk
    if 'make_bulk_list' in args.stage:
        ukb_preprocess.make_bulk_subj_list(config_json, args)

    # download bulk field
    if 'download_bulk' in args.stage:
        ukb_preprocess.download_bulk_stage(config_json, args)

    # download genetic data
    if 'download_genetic' in args.stage:
        ukb_preprocess.download_genetic(config_json, args)

    # compile already downloaded bulk data into a single csv
    if 'bulk_to_csv' in args.stage:
        ukb_preprocess.bulk_to_csv(config_json, args)

    # RUN the phesant pipeline to precompute AGE regressions
    if 'run_phesant' in args.stage:
        ukb_phesant.run_phesant(config_json, args)

    # PHESANT output comes in multiple parts, combine them into a single file
    if 'prep_data_for_phesant' in args.stage:
        ukb_phesant.prep_data_for_phesant(config_json, args)

    # create the combined UKB SQL database
    if 'ukb_sql' in args.stage:
        ukb_sql.ukb_create_sql(config_json, args)

    # create metadata file
    if 'make_metadata' in args.stage:
        ukb_sql.make_metadata(config_json, args)

    # preprocess genetic data
    if 'snp_preprocess' in args.stage:
        ukb_genetics.snp_preprocess(config_json, args)


if __name__ == '__main__':
    sys.exit(main())






