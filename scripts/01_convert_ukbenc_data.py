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
    parser.add_argument('--overwrite', dest='overwrite', default=False, required=False, action='store_true', help='overwrite')
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
            if os.path.exists(dest_path):
                os.remove(dest_path)
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
        converted_data = './{}.enc_ukb'.format(enc_base)
        if not os.path.exists(converted_data) or opt.overwrite == True:
            if os.path.exists(converted_data):
                os.remove(converted_data)
            subprocess.call([os.path.join(config_json['repo_dir'], 'external/ukbunpack'), './{}.enc'.format(enc_base), './{}'.format(key_name)])
            #subprocess.call(['/ukbtools/ukbunpack', './{}.enc'.format(enc_base), './{}'.format(key_name)])

        # unpack data - r
        tab_file = './{}.tab'.format(enc_base)
        if not os.path.exists(tab_file) or opt.overwrite == True:
            if os.path.exists(tab_file):
                os.remove(tab_file)
            subprocess.call([os.path.join(config_json['repo_dir'], 'external/ukbconv'), './{}.enc_ukb'.format(enc_base), 'r', '-o{}'.format(enc_base), '-e{}'.format('/ref_files/encoding.ukb')])
            #subprocess.call(['/ukbtools/ukbconv', './{}.enc_ukb'.format(enc_base), 'r', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])

        # unpack data - csv
        csv_file = './{}.csv'.format(enc_base)
        if not os.path.exists(csv_file) or opt.overwrite == True:
            if os.path.exists(csv_file):
                os.remove(csv_file)
            subprocess.call([os.path.join(config_json['repo_dir'], 'external/ukbconv'), './{}.enc_ukb'.format(enc_base), 'csv', '-o{}'.format(enc_base), '-e{}'.format('/ref_files/encoding.ukb')])
            #subprocess.call(['/ukbtools/ukbconv', './{}.enc_ukb'.format(enc_base), 'csv', '-ot{}'.format(enc_base), '-e{}'.format('./encoding.ukb')])



if __name__ == "__main__":
    main()


















