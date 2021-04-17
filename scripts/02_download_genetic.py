#!/bin/python

import os
import json
import subprocess
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    parser.add_argument('--download_genotyped', dest='download_genotyped', action='store_true', default=False, required=False)
    parser.add_argument('--slurm', '-e', dest='slurm', action='store_true', default=False)
    opt = parser.parse_args()

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.download_genotyped = True
        opt.slurm = True

    # list of chromosomes to download
    chrom_list = [str(chr) for chr in range(1,27)]
    chrom_list = chrom_list + ['X', 'Y', 'XY', 'MT']


    # copy encoding files/keys to the genotyped data directories
    enc_file      = config_json['ukb_enc'].replace('.enc', '.enc_ukb')
    enc_base      = enc_file.split('/')[-1]
    enc_key       = config_json['ukb_key']
    key_base      = enc_key.split('/')[-1]

    genotyped_dir = os.path.join(config_json['base_dir'], 'data/ukb/genotypes/genotyped')
    if not os.path.exists(os.path.join(genotyped_dir, key_base)):
        shutil.copy(enc_key, os.path.join(genotyped_dir, key_base))

    imputed_dir = os.path.join(config_json['base_dir'], 'data/ukb/genotypes/imputed')
    if not os.path.exists(os.path.join(imputed_dir, key_base)):
        shutil.copy(enc_key, os.path.join(imputed_dir, key_base))

    # copy the ukbgene utility
    orig_ukbgene = os.path.join(config_json['base_dir'], 'external/ukbgene')
    shutil.copy(orig_ukbgene, os.path.join(genotyped_dir, 'ukbgene'))
    shutil.copy(orig_ukbgene, os.path.join(imputed_dir, 'ukbgene'))

    if opt.download_genotyped == True:
        for chrom in chrom_list:
            os.chdir(genotyped_dir)
            os.system()
            fetch_cmd = ['./ukbgene', 'cal', '-a./{}'.format(key_base), '-c{}'.format(chrom)]
            if slurm == True:
                

            p = subprocess.Popen(fetch_cmd, stdout=subprocess.PIPE)
            out, err = p.communicate()

            print(chrom)


    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.slurm = True

    #'./ukbgene cal -a../../raw/ukb40501.key -c21'

    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    ukbgene typename -cchrom [flags]




