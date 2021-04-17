#!/bin/python

import os
import json
import shutil
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', action='config', required=True)
    opt = parser.parse_args()
    config_file = opt.config

    # read user configuration file
    #config_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # path to UKB encoding file
    ukb_enc  = config_json['ukb_enc']
    ukb_key  = config_json['ukb_key']
    enc_name = config_json['ukb_enc'].split('/')[-1]
    key_name = config_json['ukb_key'].split('/')[-1]

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
    enc_dir  = '/'.join(ukb_enc.split('/')[:-1])
    enc_name = ukb_enc.split('/')[-1]
    enc_base = enc_name.split('.')[0]

    # check to see if the ukbconv utility is downloaded
    #ukb_conv_path = os.path.join(enc_dir, 'ukbconv')
    #if not os.path.exists(ukb_conv_path):
    #    subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/ukbconv', '--directory-prefix={}'.format(enc_dir)])

    # download encoding file, if needed
    #encoding_file = os.path.join(enc_dir, 'encoding.ukb')
    #if not os.path.exists(encoding_file):
    #    subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb', '--directory-prefix={}'.format(enc_dir)])

    # cd to data directory
    os.chdir(enc_dir)

    # convert command
    subprocess.call(['/app/ukbconv', enc_name, 'r', '-otest{}'.format(enc_base), '-e/ref_files/encoding.ukb'])


if __name__ == "__main__":
    main()


















