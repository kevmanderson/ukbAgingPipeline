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
import zipfile
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
    log.info('RUNNING: STAGE = RUN_PHESANT')

    parser = argparse.ArgumentParser()
    parser.add_argument('--sumstat_in', dest='sumstat_in', required=True, default=None)
    parser.add_argument('--sumstat_out', dest='sumstat_out', required=True, default=None)
    parser.add_argument('--ldsc_munge_py', dest='ldsc_munge_py', required=True, default=None)

    parser = argparse.ArgumentParser()
    args   = parser.parse_args()
    args.sumstat_in    = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw/IGAP_stage_1.txt'
    args.out_dir   = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/munged/'
    args.out_name  = 'Alz_test'
    args.ldsc_munge_py = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/ldsc/munge_sumstats.py'
    args.a1 = 'Non_Effect_Allele'
    args.a2 = 'Effect_allele'
    args.p = 'Pvalue'
    args.frq = None
    args.n_cas = 17008
    args.n_con = 37154
    args.n_cas_col = 'Ncases'
    args.n_con_col = 'Ncontrols'
    args.daner_n = False
    args.signed_sumstats = 'Beta'
    args.snp = 'MarkerName'
    
    out_path = os.path.join(args.out_dir, args.out_name)

    # LDSC munge_stats utility
    munge_base = f'''
source activate ldsc

cd {args.out_dir}

{args.ldsc_munge_py} \\
    --sumstats={args.sumstat_in} \\
    --out={out_path} \\
    --a1={args.a1} \\
    --a2={args.a2} \\
    --p={args.p} \\
    --signed-sumstats={args.signed_sumstats},0 \\
    --snp={args.snp} \\'''

    # add N column
    if args.n_cas != None and args.n_con != None:
        munge_base = f'''{munge_base}
--N-cas={args.n_cas} \\
--N-con={args.n_con}
        '''
    elif args.n_cas_col != None and args.n_con_col != None:
        munge_base = f'''{munge_base}
--N-cas-col={args.n_cas_col} \\
--N-con-col={args.n_con_col}
        '''
    print(munge_base)

    file_path = args.sumstat_in
    unzip_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/munged'


def unzip_sumstat(file_path, unzip_dir):
    if '.zip' in file_path.split('/')[-1]:
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall(unzip_dir)




if __name__ == "__main__":
    main()

















