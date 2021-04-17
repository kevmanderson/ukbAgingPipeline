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
    # ------
    #config_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    base_dir = config_json['base_dir']

    ref_dir = os.path.join(base_dir, 'ref_files')
    if not os.path.exists(ref_dir):
        os.mkdir(ref_dir)

    data_dir = os.path.join(base_dir, 'data')
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    ukb_dir = os.path.join(data_dir, 'ukb')
    if not os.path.exists(ukb_dir):
        os.mkdir(ukb_dir)
    raw_dir = os.path.join(data_dir, 'raw')
    if not os.path.exists(raw_dir):
        os.mkdir(raw_dir)
    bulk_dir = os.path.join(data_dir, 'mri_bulk')
    if not os.path.exists(bulk_dir):
        os.mkdir(bulk_dir)
    gene_dir = os.path.join(data_dir, 'genetics')
    if not os.path.exists(gene_dir):
        os.mkdir(gene_dir)

if __name__ == "__main__":
    main()



