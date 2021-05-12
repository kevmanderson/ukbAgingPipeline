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
import argparse
import numpy as np

from utilities.utilities import create_directories, write_slurm, submit_slurm
from download import decrypt_ukb_data, get_ukbutils


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
        args.stage = 'prep'
        args.slurm_partition = 'short'
        args.stages = ['decrypt']
        args.convert_all_fields = false
        args.singularity_container = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/simons_ukb_aging_pipeline'
        return args

    else:
        parser = argparse.ArgumentParser()
        parser.add_argument('--config', '-c', dest='config', required=True)
        parser.add_argument('--stage', dest='stage', required=False)
        parser.add_argument('--bulk-field', '-b', dest='bulk_field', action='append', nargs='+', required=False)
        parser.add_argument('--make-bulk-list', dest='make_bulk_list', action='store_true', default=False)
        parser.add_argument('--download-bulk-data', dest='download_bulk_data', action='store_true', default=False)
        parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
        parser.add_argument('--slurm_partition', '-p', dest='slurm_partition', required=False, default=False)
        parser.add_argument('--singularity_container', dest='singularity_container', required=False, default=False)
        parser.add_argument('--stages', dest='stages', nargs='+', required=False, default=False)
        parser.add_argument('--convert_all_fields', action='store_true', required=False, default=False)

        args = parser.parse_args()
        return args


def main(argv=None):

    args = create_parser()
    #args = create_parser('yale')
    print(args)

    # read configuration file
    with open(args.config, 'r') as f:
        config_json = json.load(f)[0]
    slurm_dir = os.path.join(config_json['base_dir'], 'slurm')

    # create project directory structure if necessary
    create_directories(root_dir=config_json['base_dir'])

    # download UKB utilities if needed
    get_ukbutils(util_dir=os.path.join(config_json['base_dir'], 'data/ukb/raw'))

    # decrypt the encoded ukb data
    if 'decrypt' in args.stages:
        decrypt_stage(config_json, args)

    # convert encoded data files
    if 'convert' in args.stages:
        convert_stage(config_json, args)


def convert_stage(config_json, args):
    '''
    Convert decrypted UKB data to readable formats
    '''

    # ukb utility to convert encrypted data
    ukbconv = os.path.join(config_json['base_dir'], 'data/ukb/external/ukbconv')

    if args.convert_all_fields == False:

        cat_file  = os.path.join(config_json['repo_dir'], 'ref_files/funpack/categories.tsv')
        field_arr = read_funpack_categories(cat_file)

        # write list of fields to convert
        raw_dir    = os.path.join(config_json['base_dir'], 'data/ukb/raw')
        field_file = os.path.join(raw_dir, 'fields_to_convert.txt')
        field_arr  = list(np.sort(field_arr))
        with open(field_file, 'w') as f:
            [f.writelines('{}\n'.format(field) for field in field_arr)]

        # decrypt all input enc/key pairs
        for ukb_enc in config_json['ukb_encs']:
            enc_ukb_file = ukb_enc['ukb_enc'].split('/')[-1].replace('.enc','.enc_ukb')
            conv_cmd = convert_cmd(ukbconv=ukbconv, enc_ukb=enc_ukb_file)

            decrypt_cmd = f"cd {raw_dir}\n"
            decrypt_cmd = f"{decrypt_cmd}{ukbunpack} ./{enc_name} ./{key_name}"


def convert_cmd(ukbconv, raw_dir, enc_ukb, output_format):
    '''Create command to convert ukb data fields'''

    decrypt_cmd = f"cd {raw_dir}\n"
    decrypt_cmd = f"{decrypt_cmd}{ukbconv} {enc_ukb} {output_format}./{enc_name} ./{key_name}"


def decrypt_stage(config_json, args):
    '''
    Decrypt encoded ukb data stage
    '''

    # path where data will be symlinked/decrypted
    raw_dir   = os.path.join(config_json['base_dir'], 'data/ukb/raw')
    write_dir = os.path.join(config_json['base_dir'], 'slurm')

    # ukb utility to unpack encrypted data
    ukbunpack = os.path.join(config_json['base_dir'], 'data/ukb/external/ukbunpack')

    # decrypt all input enc/key pairs
    for ukb_enc in config_json['ukb_encs']:
        enc_base = ukb_enc['ukb_enc'].split('/')[-1].replace('.enc', '')

        # create the ukb file decryption command
        decrypt_cmd = decrypt_ukb_data(ukbunpack=ukbunpack, write_dir=raw_dir, enc=ukb_enc['ukb_enc'], key=ukb_enc['ukb_key'])

        write_file = os.path.join(write_dir, 'decrypt_{}'.format(enc_base))
        if args.slurm == True:
            submit_slurm(write_slurm(slurm_file=write_file, partition=args.slurm_partition, cmd=decrypt_cmd, jobName=enc_base, stime='6:00:00', nthreads=2))
        else:
            write_cmd(cmd=decrypt_cmd, write_file='{}.bash'.format(write_file), print_path=True)



if __name__ == '__main__':
    sys.exit(main())






