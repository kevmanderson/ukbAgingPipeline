

import os
import json
import subprocess
import argparse
import glob

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

    print('Configuration file:')
    print(print(json.dumps(config_json, indent=2)))

    base_dir = config_json["base_dir"]
    print('Base Directory: {}'.format(base_dir))

    data_dir = os.path.join(base_dir, 'data')
    print('Creating Data Directory: {}'.format(data_dir))
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    ukb_dir = os.path.join(data_dir, 'ukb')
    print('Creating UKB Directory: {}'.format(ukb_dir))
    if not os.path.exists(ukb_dir):
        os.mkdir(ukb_dir)

    bulk_dir = os.path.join(ukb_dir, 'bulk')
    print('Creating Bulk Directory: {}'.format(bulk_dir))
    if not os.path.exists(bulk_dir):
        os.mkdir(bulk_dir)

    raw_dir = os.path.join(ukb_dir, 'raw')
    print('Creating Raw Directory: {}'.format(raw_dir))
    if not os.path.exists(raw_dir):
        os.mkdir(raw_dir)

    slurm_dir = os.path.join(base_dir, 'slurm')
    print('Creating Slurm Directory: {}'.format(slurm_dir))
    if not os.path.exists(slurm_dir):
        os.mkdir(slurm_dir)

    gene_dir = os.path.join(ukb_dir, 'genotypes')
    print('Creating Genetics Directory: {}'.format(gene_dir))
    if not os.path.exists(gene_dir):
        os.mkdir(gene_dir)

    genotyped_dir = os.path.join(gene_dir, 'genotyped')
    print('Creating Genetics Directory: {}'.format(genotyped_dir))
    if not os.path.exists(genotyped_dir):
        os.mkdir(genotyped_dir)

    imputed_dir = os.path.join(gene_dir, 'imputed')
    print('Creating Genetics Directory: {}'.format(imputed_dir))
    if not os.path.exists(imputed_dir):
        os.mkdir(imputed_dir)

    external_dir = os.path.join(data_dir, 'external')
    print('Creating External Directory: {}'.format(external_dir))
    if not os.path.exists(external_dir):
        os.mkdir(external_dir)
    for filename in glob.glob(os.path.join('/opt/ukbtools', 'ukb*')):
        shutil.copy(filename, external_dir)



if __name__ == "__main__":
    main()

    