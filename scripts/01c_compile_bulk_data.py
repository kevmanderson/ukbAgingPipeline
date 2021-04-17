#!/bin/python

import os
import json
import subprocess
import argparse
import pandas as pd
import glob
import numpy as np
from functools import reduce
import sys
sys.path.append('/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts')
sys.path.append('/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/utilities')
from utilities import submitSlurm, writeSlurm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    parser.add_argument('--bulk-field', '-b', dest='bulk_field', action='append', nargs='+' required=True)
    parser.add_argument('--slurm', '-e', dest='slurm', action='store_true', default=False)
    opt = parser.parse_args()

    print(opt.bulk_field)

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.bulk_id = 25755
        opt.num_components = 100
        opt.is_rsfa = True
        opt.is_funconn = False
        opt.slurm = True

    config_file = opt.config
    bulk_id     = opt.bulk_id
    slurm       = opt.slurm
    num_components = opt.num_components
    is_funconn = opt.is_funconn
    is_rsfa = opt.is_rsfa

    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]


    component_path = '/ref_files/rfMRI_GoodComponents_d{}_v1.txt'.format(num_components)
    component_in   = pd.read_csv(component_path, delim_whitespace=True, header=None)
    component_map  = {i: x for i, x in enumerate(component_in.transpose()[0])}

    # directory where all the operations will take place
    ukb_enc = config_json['ukb_enc']
    enc_dir = '/'.join(ukb_enc.split('/')[:-1])

    # create directory for all the output bulk files
    bulk_dir   = os.path.join(enc_dir, 'bulk_{}'.format(bulk_id))
    bulk_files = glob.glob(os.path.join(bulk_dir, '*{}*.txt'.format(bulk_id)))

    bulk_file = bulk_files[0]
    is_rsfa = True

    if is_rsfa:
        bulk_dict = {}
        for i,bulk_file in enumerate(bulk_files):
            if i % 100 == 0:
                print(i)
            visit = bulk_file.split('/')[-1].split('_')[-2]
            subj  = bulk_file.split('/')[-1].split('_')[0]
            if visit not in bulk_dict.keys():
                bulk_dict[visit] = []
            bulk_data  = pd.read_csv(bulk_file, delim_whitespace=True, header=None)
            bulk_names = ['{}_{}_{}'.format(bulk_id, visit, component_map[i]) for i in range(bulk_data.shape[1])]
            bulk_data.columns = bulk_names
            bulk_data.insert(0, 'id', subj)
            bulk_dict[visit].append(bulk_data)
        df_list = [pd.concat(bulk_dict[key]) for key in bulk_dict.keys()]

    df_write = os.path.join(bulk_dir, )
    df = reduce(lambda df1,df2: pd.merge(df1, df2, on='id', how='outer'), df_list)

if __name__ == "__main__":
    main()


