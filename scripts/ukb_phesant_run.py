#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# Convert decrypted UKB data to readable formats

import os
import glob
import datetime
import calendar
import argparse
import logging
import numpy as np
import pandas as pd
import textwrap
from functools import reduce
from argparse import RawTextHelpFormatter
import utilities.utilities as utilities
log = logging.getLogger('ukb')

def main():
    '''
    Before downloading bulk, data first need to make lists of subjecst that possess data

    Parameters
    ----------


    Returns
    -------

    '''

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    log.info('RUNNING: STAGE = RUN_PHESANT')

    parser = argparse.ArgumentParser()
    parser.add_argument('--pheno_file', dest='pheno_file', required=True, default=None,
                        help=textwrap.dedent('''
                         PHESANT formatted phenotype file
                        '''))
    parser.add_argument('--confounder_file', dest='confounder_file', required=True, default=None,
                        help=textwrap.dedent('''
                         PHESANT formatted confounder file
                        '''))
    parser.add_argument('--age_var', dest='age_var', required=True, default=None,
                        help='name of the age variable')
    parser.add_argument('--visit', dest='visit', required=True, default=None,
                        help='name of the age variable')
    parser.add_argument('--datacoding_file', dest='datacoding_file', required=True, default=None,
                        help='Ordinal data codings')
    parser.add_argument('--out_dir', dest='out_dir', required=True, default=None,
                        help='output directory')
    parser.add_argument('--nparts', dest='nparts', required=True, default=None, type=int,
                        help='nparts')
    parser.add_argument('--metadata_file', dest='metadata_file', required=True, default=None,
                        help='metadata_file')
    parser.add_argument('--slurm', dest='slurm', required=False, default=False,
                        help='submit to slurm or not')
    parser.add_argument('--write_scripts_dir', dest='write_scripts_dir', required=False, default=False,
                        help='submit to slurm or not')
    parser.add_argument('--phesant_dir', dest='phesant_dir', required=False, default=False,
                        help='submit to slurm or not')
    args = parser.parse_args()
    print(args)

    # set up filepaths
    #phesant_dir = os.path.join(config_json['repo_dir'], 'external/PHESANT/WAS')
    #resDir = os.path.join(config_json['base_dir'], 'data/ukb/phesant/')
    #variablelistfile = os.path.join(config_json['repo_dir'], 'ref_files/aging-outcome-info.tsv')
    #datacodingfile = os.path.join(config_json['repo_dir'], 'ref_files/data-coding-ordinal-info.txt')

    #nparts = 30
    #for visit in visit_arr:
        #age_var = 'x21003_{}_0'.format(visit)
        #confounderfile = os.path.join(config_json['base_dir'], 'data/ukb/raw',
        #                                'phesant_covar_vis{}.csv'.format(visit))
        #phenofile = os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_visit{}.csv'.format(visit))

    for part in range(1, args.nparts + 1):
        phesant_cmd = f'''#!/bin/bash
cd {args.phesant_dir}
Rscript phenomeScan.r \\
        --phenofile={args.pheno_file} \\
        --variablelistfile={args.metadata_file} \\
        --confounderfile={args.confounder_file} \\
        --datacodingfile={args.datacoding_file} \\
        --traitofinterest="{args.age_var}" \\
        --resDir={args.out_dir} \\
        --userId="xeid" \\
        --sensitivity \\
        --numParts={args.nparts} --partIdx={part} \\
        --ageistraitofinterest \\
        --genetic="FALSE" \\
        --standardise="TRUE" \\
        --visit={args.visit} \\
        --file_prepend=phewas_visit{args.visit} \\
        --save'''

        # write executable to file
        if args.slurm == True:
            log.info('submit slurm')
            write_file = os.path.join(args.write_scripts_dir, 'phesant_visit{}_part{}-{}'.format(args.visit, part, args.nparts))
            utilities.submit_slurm(
                utilities.write_slurm(slurm_file=write_file,
                                        partition=args.slurm_partition, cmd=phesant_cmd,
                                        jobName='{}_{}'.format(visit, part), stime='6:00:00', nthreads=12))
        else:
            write_file = os.path.join(args.write_scripts_dir, 'phesant_visit{}_part{}-{}.bash'.format(args.visit, part, args.nparts))
            utilities.write_cmd(cmd=phesant_cmd, write_file=write_file)


if __name__ == "__main__":
    main()
