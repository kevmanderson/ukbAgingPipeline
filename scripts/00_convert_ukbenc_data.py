#!/bin/python

import os
import json
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    opt = parser.parse_args()
    config_file = opt.config

    # read user configuration file
    #config_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # path to UKB encoding file
    ukb_enc = config_json['ukb_enc']

    # directory where all the operations will take place
    enc_dir = '/'.join(ukb_enc.split('/')[:-1])
    enc_name = ukb_enc.split('/')[-1]
    enc_base = enc_name.split('.')[0]

    # check to see if the ukbconv utility is downloaded
    ukb_conv_path = os.path.join(enc_dir, 'ukbconv')
    if not os.path.exists(ukb_conv_path):
        subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/ukbconv', '--directory-prefix={}'.format(enc_dir)])

    # download encoding file, if needed
    encoding_file = os.path.join(enc_dir, 'encoding.ukb')
    if not os.path.exists(encoding_file):
        subprocess.call(['wget', '-nd', 'biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb', '--directory-prefix={}'.format(enc_dir)])

    # cd to data directory
    os.chdir(enc_dir)

    # convert r command
    subprocess.call([ukb_conv_path, enc_name, 'r', '-o{}'.format(enc_base), '-eencoding.ukb'])

    # convert csv command (for PHESANT)
    subprocess.call([ukb_conv_path, enc_name, 'csv', '-o{}'.format(enc_base), '-eencoding.ukb'])


if __name__ == "__main__":
    main()



