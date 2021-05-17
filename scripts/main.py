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

import utilities.utilities as utilities

from utilities.utilities import write_slurm, submit_slurm, make_dir, make_symlink
from download import decrypt_ukb_data, get_ukbutils


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
        parser.add_argument('--log_file', dest='log_file', default=None)
        parser.add_argument('--quiet', '-q', action='store_true', default=False)
        parser.add_argument('--noisy', '-n', action='count', default=0)

        args = parser.parse_args()
        return args


def main(argv=None):

    logging.basicConfig(level=logging.DEBUG)

    date = datetime.date.today()

    args = create_parser()
    # args = create_parser('yale')
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
        decrypt_stage(config_json, args)

    # convert encoded data files
    if 'convert' in args.stage:
        convert_stage(config_json, args)

    # make list of subjects to download for bulk
    if 'make_bulk_list' in args.stage:
        make_bulk_subj_list(config_json, args)

    # download bulk field
    if 'download_bulk' in args.stage:
        download_bulk_stage(config_json, args)

    # compile already downloaded bulk data to csv
    if 'bulk_to_csv' in args.stage:
        bulk_to_csv(config_json, args)



def bulk_to_csv(config_json, args):
    '''Concat all bulk data into a single csv'''
    log.info('RUNNING: STAGE = BULK_TO_CSV')

    # point to directories
    bulk_dir  = os.path.join(config_json['base_dir'], 'data/ukb/bulk')
    raw_dir   = os.path.join(config_json['base_dir'], 'data/ukb/raw')

    # iterate through all of the bulk fields to download
    for field in args.bulk_field:
        # get the descriptive bulk_name and ukb bulk_id
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]

        # where the bulk data live
        bulk_field_dir = os.path.join(bulk_dir, bulk_name)
        log.debug('Read BULK_ID ("{}") from DIR ("{}")'.format(bulk_id, bulk_field_dir))

        # get paths to downloaded files
        bulk_path_df = {}
        bulk_path_df['2'] = glob.glob(os.path.join(bulk_field_dir, '*{}*_2_0.txt'.format(bulk_id)))
        bulk_path_df['3'] = glob.glob(os.path.join(bulk_field_dir, '*{}*_3_0.txt'.format(bulk_id)))

        # 25 component features
        if bulk_id in ['25750', '25752', '25754']:
            n_comp = 21
            with open(os.path.join(config_json['repo_dir'], 'ref_files', 'rfMRI_GoodComponents_d25_v1.txt'), 'r') as f: comp_names = f.readline().split()
        # 100 component features
        elif bulk_id in ['25751', '25753', '25755']:
            n_comp = 55
            with open(os.path.join(config_json['repo_dir'], 'ref_files', 'rfMRI_GoodComponents_d100_v1.txt'), 'r') as f: comp_names = f.readline().split()


        #if bulk_id in ['25754', '25755']:
        for vis in ['2', '3']:
            # RSFA
            if bulk_id in ['25754', '25755']:
                n_values  = n_comp
                col_names = ['{}_{}_r{}'.format(bulk_name, vis, x) for x in comp_names]
            # func-con
            elif bulk_id in ['25750', '25751', '25752', '25753']:
                n_values = ((n_comp * n_comp) - n_comp) / 2
                zero_df = pd.DataFrame(np.zeros((n_comp, n_comp)))
                zero_df.columns = comp_names
                zero_df.index   = comp_names
                zero_df = zero_df.where(np.triu(np.ones(zero_df.shape), k=1).astype(np.bool))
                zero_df = zero_df.stack().reset_index()
                edges = zip(zero_df[zero_df.columns[0]], zero_df[zero_df.columns[1]])
                col_names = ['{}_{}_r{}r{}'.format(bulk_name, vis, x[0], x[1]) for x in edges]

            # read the individual data files
            data_list  = []
            bulk_paths = bulk_path_df[vis]
            for i,bulk_path in enumerate(bulk_paths):
                id = bulk_path.split('/')[-1].split('_')[0]
                #if i % 10000 == 0: print(i)
                with open(bulk_path) as f:
                    arr = f.readline()
                data_row = np.array(arr.split()).astype(float).astype(str)
                if len(data_row) == n_values:
                    keep_row = [id] + list(data_row)
                    data_list.append(keep_row)

            # organize into dataframe
            bulk_df = pd.DataFrame(data_list)
            bulk_df.columns = ['id'] + col_names

            # write bulk data to csv
            csv_path = os.path.join(raw_dir, 'bulk_{}_{}.csv'.format(bulk_id, vis))
            log.info('Writing bulk data to: {}'.format(csv_path))
            bulk_df.to_csv(csv_path, index=None)


def make_bulk_subj_list(config_json, args):
    '''Before downloading bulk, data first need to make lists of subjecst that possess data'''
    log.info('RUNNING: STAGE = MAKE_BULK_LIST')

    ukbconv = os.path.join(config_json['base_dir'], 'data/ukb/external/ukbconv')

    # point to directories
    raw_dir   = os.path.join(config_json['base_dir'], 'data/ukb/raw')
    bulk_dir  = os.path.join(config_json['base_dir'], 'data/ukb/bulk')
    write_dir = os.path.join(config_json['base_dir'], 'slurm')

    # iterate through all of the bulk fields to download
    for field in args.bulk_field:
        # get the descriptive bulk_name and ukb bulk_id
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]

        # create bulk output directory if it doesnt exist
        bulk_out_dir = os.path.join(bulk_dir, bulk_name)
        make_dir(bulk_out_dir)

        # there are often more than one data bucket, so check to see which one contains the bulk id of interest
        for ukb_enc_l in config_json['ukb_encs']:
            # path to the converted ukb csv file
            enc_base     = ukb_enc_l['ukb_enc'].split('/')[-1].replace('.enc', '')
            ukb_enc      = os.path.join(raw_dir, '{}.csv'.format(enc_base))
            ukb_hdr_cols = utilities.get_ukbenc_header(ukb_enc)

            # if bulk id is in this basket
            if bulk_id in ukb_hdr_cols:
                log.debug('BULK ID ("{}") is in BUCKET ("{}")'.format(bulk_id, enc_base))
                enc_ukb_use  = '{}.enc_ukb'.format(enc_base)
                ukb_enc_src  = os.path.join(raw_dir, enc_ukb_use)
                ukb_enc_dest = os.path.join(bulk_out_dir, enc_ukb_use)
                make_symlink(ukb_enc_src, ukb_enc_dest)

        # build command for slurm submission
        convert_cmd  = 'cd {}'.format(bulk_out_dir)
        download_cmd = [ukbconv, './{}'.format(enc_ukb_use), 'bulk', '-s{}'.format(bulk_id), '-o./{}'.format(bulk_id)]
        submit_cmd   = '{}\n{}'.format(convert_cmd, ' '.join(download_cmd))

        write_file = os.path.join(write_dir, 'make_bulk_sublist_{}.bash'.format(bulk_id))
        utilities.write_cmd(cmd=submit_cmd, write_file=write_file)




def download_bulk_stage(config_json, args):
    '''Create scripts to download bulk data'''

    ukbfetch  = os.path.join(config_json['base_dir'], 'data/ukb/external/ukbfetch')
    bulk_dir  = os.path.join(config_json['base_dir'], 'data/ukb/bulk')
    write_dir = os.path.join(config_json['base_dir'], 'slurm')

    print('Creating Bulk Files for Download')
    for field in args.bulk_field:
        # get the descriptive bulk_name and ukb bulk_id
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]

        # where the bulk data live
        bulk_out_dir = os.path.join(bulk_dir, bulk_name)

        # TODO: hacky but ok for now
        key_path = glob.glob(os.path.join(bulk_out_dir, '*key'))[0]
        key_file = key_path.split('/')[-1]

        # copy the ukbfetch utility to the bulk download dir
        make_symlink(ukbfetch, os.path.join(bulk_out_dir, 'ukbfetch'))

        # read file with subjects to pull
        bulk_file = os.path.join(bulk_out_dir, '{}.bulk'.format(bulk_id))
        bulk_df = pd.read_csv(bulk_file, header=None, delim_whitespace=True)
        downloaded_paths  = glob.glob(os.path.join(bulk_out_dir, '*{}*'.format(bulk_id)))
        downloaded_files  = set([x.split('/')[-1] for x in downloaded_paths])
        total_files       = ['{}_{}.txt'.format(i[0], i[1]) for i in zip(bulk_df[0], bulk_df[1])]
        files_to_download = list(set(total_files) - set(downloaded_files))

        log.debug('Reading bulk subjects: {}'.format(bulk_file))
        log.debug('{} total files'.format(len(total_files)))
        log.debug('{} downloaded files'.format(len(downloaded_files)))
        log.debug('{} files to download'.format(len(files_to_download)))

        if len(files_to_download) > 0:
            download_bulk_file = bulk_file.replace('.bulk', '_use.bulk')
            subs_to_download   = list(set([x.split('_')[0] for x in files_to_download]))
            bulk_download_df   = bulk_df.loc[bulk_df[0].astype(str).isin(subs_to_download)]
            bulk_download_df.to_csv(download_bulk_file, sep='\t', header=None, index=None)

            start = 1
            nrows = bulk_download_df.shape[0]
            job_array = []
            while start < nrows:
                print(start)
                slurm_file = os.path.join(write_dir, 'download_bulk_id{}_s{}.bash'.format(bulk_id, start))
                job_array.append(slurm_file)

                # fetch command
                fetch_cmd = f'''cd {bulk_out_dir}'''
                fetch_cmd = f'''{fetch_cmd}\n\n./ukbfetch \\\n {'-b./{}_use.bulk'.format(bulk_id)} \\\n -a./{key_file} \\\n -s{start} -m5000'''
                start = start + 5000

                if args.slurm == True:
                    slurm_path = writeSlurm(slurm_file, slurm_partition, fetch_cmd, str(start), stime='6:00:00', n_gpu=None,
                                            nthreads=2, mem='8G')
                    print(slurm_path)
                    job_id = submitSlurm(slurm_path)

                else:
                    utilities.write_cmd(cmd=fetch_cmd, write_file=slurm_file)


                    # write script to file
                    with open(slurm_file, 'w') as f:
                        f.write(fetch_cmd)
                    os.system('chmod 770 {}'.format(slurm_file))


def convert_stage(config_json, args):
    '''Convert decrypted UKB data to readable formats'''

    # ukb utility to convert encrypted data
    log.info('RUNNING: STAGE = CONVERT')
    ukbconv   = os.path.join(config_json['base_dir'], 'data/ukb/external/ukbconv')
    write_dir = os.path.join(config_json['base_dir'], 'slurm')

    if args.convert_all_fields == False:

        cat_file  = os.path.join(config_json['repo_dir'], 'ref_files/funpack/categories.tsv')
        field_arr = utilities.read_funpack_categories(cat_file)

        # write list of fields to convert
        raw_dir    = os.path.join(config_json['base_dir'], 'data/ukb/raw')
        field_file = os.path.join(raw_dir, 'fields_to_convert.txt')
        field_arr  = list(np.sort(list(set(field_arr))))
        with open(field_file, 'w') as f:
            [f.writelines('{}\n'.format(field) for field in field_arr)]

        # decrypt all input enc/key pairs
        output_format = 'csv'
        for ukb_enc in config_json['ukb_encs']:
            for output_format in ['csv', 'r']:
                ukb_key      = ukb_enc['ukb_key'].split('/')[-1]
                enc_ukb_file = ukb_enc['ukb_enc'].split('/')[-1].replace('.enc','.enc_ukb')
                conv_cmd = convert_cmd(ukbconv=ukbconv, raw_dir=raw_dir, enc_ukb=enc_ukb_file, output_format=output_format, field_file='fields_to_convert.txt')

                # write command
                write_file = os.path.join(write_dir, 'convert_{}_{}.bash'.format(enc_ukb_file.replace('.enc_ukb',''), output_format))
                utilities.write_cmd(cmd=conv_cmd, write_file=write_file)


def convert_cmd(ukbconv, raw_dir, enc_ukb, output_format, field_file=None):
    '''Create command to convert ukb data fields'''
    o_name = enc_ukb.replace('.enc_ukb','')
    decrypt_cmd = f"cd {raw_dir}\n"
    if field_file != None:
        decrypt_cmd = f"{decrypt_cmd}{ukbconv} {enc_ukb} {output_format} -o./{o_name} -i./{field_file}"
    else:
        decrypt_cmd = f"{decrypt_cmd}{ukbconv} {enc_ukb} {output_format} -o./{o_name}"
    return decrypt_cmd


def decrypt_stage(config_json, args):
    '''
    Decrypt encoded ukb data stage
    '''
    log.info('RUNNING: STAGE = DECRYPT')

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
            utilities.submit_slurm(utilities.write_slurm(slurm_file=write_file, partition=args.slurm_partition, cmd=decrypt_cmd, jobName=enc_base, stime='6:00:00', nthreads=2))
        else:
            utilities.write_cmd(cmd=decrypt_cmd, write_file='{}.bash'.format(write_file))



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






