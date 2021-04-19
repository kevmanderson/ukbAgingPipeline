#!/bin/python

import os
import json
import shutil
import subprocess
import argparse

def main():
    example_text = '''example: 
     python test.py -t template/test.py'''

    parser = argparse.ArgumentParser(epilog=example_text, add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', '-c', dest='config', required=True, help='Full path to the user configuration file')
    parser.add_argument('--r_only', dest='r_only', default=False, required=False, action='store_true', help='Only convert data into the r-readable format')
    parser.add_argument('--csv_only', dest='csv_only', default=False, required=False, action='store_true', help='Only convert data into csv format')
    parser.add_argument('--help', '-h', action='help', default=argparse.SUPPRESS, help='This script will unpack and convert UK Biobank')
    opt = parser.parse_args()
    config_file = opt.config

    # read user configuration file
    # config_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
    # config_file = '/ncf/sba01/ukbAgingPipeline/config.json'
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # path to UKB encoding file
    for data_pair in config_json['ukb_encs']:
        print(data_pair)
        ukb_enc  = data_pair['ukb_enc']
        ukb_key  = data_pair['ukb_key']
        enc_name = ukb_enc.split('/')[-1]
        enc_base = enc_name.replace('.enc', '')
        key_name = ukb_key.split('/')[-1]

        # path where to write all the data
        raw_dir   = os.path.join(config_json['base_dir'], 'data/ukb/raw')

        # symbolic link the data if needed
        dest_path = os.path.join(raw_dir, enc_name)
        if ukb_enc != dest_path:
            print('Symbol: {} > {}'.format(ukb_enc, dest_path))
            os.symlink(ukb_enc, dest_path)

        # copy key
        dest_key = os.path.join(raw_dir, key_name)
        if ukb_key != dest_key:
            print('Symbol: {} > {}'.format(ukb_key, dest_key))
            shutil.copyfile(ukb_key, dest_key)

        # directory where all the operations will take place
        #enc_dir  = '/'.join(ukb_enc.split('/')[:-1])
        #enc_name = ukb_enc.split('/')[-1]
        #enc_base = enc_name.split('.')[0]

        # check to see if the ukbconv utility is downloaded
        #ukb_conv_path = os.path.join(enc_dir, 'ukbconv')
        #if not os.path.exists(ukb_conv_path):
        #    subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/ukbconv', '--directory-prefix={}'.format(enc_dir)])

        # download encoding file, if needed
        #encoding_file = os.path.join(enc_dir, 'encoding.ukb')
        #if not os.path.exists(encoding_file):
        #    subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb', '--directory-prefix={}'.format(enc_dir)])

        # convert command
        print('Unpacking UKB Data')
        os.chdir(raw_dir)
        #subprocess.call(['./ukbunpack', './{}.enc'.format(enc_base), './{}'.format(key_name)])
        subprocess.call(['/app/ukbunpack', './{}.enc'.format(enc_base), './{}'.format(key_name)])

        # unpack data - r
        subprocess.call(['./ukbconv', './{}.enc_ukb'.format(enc_base), 'r', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])
        subprocess.call(['/app/ukbconv', './{}.enc_ukb'.format(enc_base), 'r', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])

        # unpack data - csv
        subprocess.call(['./ukbconv', './{}.enc_ukb'.format(enc_base), 'csv', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])
        subprocess.call(['/app/ukbconv', './{}.enc_ukb'.format(enc_base), 'csv', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])



if __name__ == "__main__":
    main()


















