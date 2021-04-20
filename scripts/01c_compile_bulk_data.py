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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    parser.add_argument('--bulk-field', '-b', dest='bulk_field', action='append', nargs='+', required=True)
    opt = parser.parse_args()

    print(opt.bulk_field)

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json'
        opt.bulk_field = [['rfmri_full_25:25750'], ['rfmri_full_100:25751'], ['rfmri_part_25:25752'], ['rfmri_part_100:25753'], ['rfmri_rsfa_25:25754'], ['rfmri_rsfa_100:25755']]

    mri_info_dict = {
        '25750': {'type': 'funcconn', 'components': 25},
        '25751': {'type': 'funcconn', 'components': 100},
        '25752': {'type': 'funcconn', 'components': 25},
        '25753': {'type': 'funcconn', 'components': 100},
        '25754': {'type': 'rsfa', 'components': 25},
        '25755': {'type': 'rsfa', 'components': 100}
    }

    config_file = opt.config
    bulk_field  = opt.bulk_field
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # dir where all the bulk data should have been written
    bulk_dir = os.path.join(config_json['base_dir'], 'data/ukb/bulk')

    for field in bulk_field:
        data_list = []
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]

        # directory and paths with the individual bulk data
        cur_bulk_dir = os.path.join(bulk_dir, bulk_name)
        bulk_files   = glob.glob(os.path.join(cur_bulk_dir, '*{}*.txt'.format(bulk_id)))

        print('ID:        {}'.format(bulk_id))
        print('Name:      {}'.format(bulk_name))
        print('Searching: {}'.format(cur_bulk_dir))
        print('Found:     {} files'.format(len(bulk_files)))

        # Resting State Functional Amplitude
        # -------------
        if mri_info_dict[bulk_id]['type'] == 'rsfa':

            # 25 or 100 spatial components
            num_components = mri_info_dict[bulk_id]['components']
            component_path = '/ref_files/rfMRI_GoodComponents_d{}_v1.txt'.format(num_components)
            component_in   = pd.read_csv(component_path, delim_whitespace=True, header=None)
            component_map  = {i: x for i, x in enumerate(component_in.transpose()[0])}

            sub_bulk_file = bulk_files[0]
            for sub_bulk_file in bulk_files:
                id = sub_bulk_file.split('/')[-1].split('_')[0]
                visit  = sub_bulk_file.split('/')[-1].split('_')[2]
                dat_in = np.loadtxt(sub_bulk_file)
                subj_data = pd.DataFrame(dat_in).transpose()
                edge_list = ['{}_{}_p{}'.format(bulk_id, visit, p) for p in component_map.values()]
                subj_data.columns = edge_list

                subj_data.insert(0, 'id', id)
                data_list.append(subj_data)

        # Resting-state Functional Connectivity
        # -------------
        if mri_info_dict[bulk_id]['type'] == 'funcconn':

            # read the mapping between parcel indices and parcel numbers (confusing)
            num_components = mri_info_dict[bulk_id]['components']
            component_path = '/ref_files/rfMRI_GoodComponents_d{}_v1.txt'.format(num_components)
            component_in   = np.loadtxt(component_path)  # pd.read_csv(component_path, delim_whitespace=True, header=None)
            component_in   = component_in.astype(int)
            parcel_arr     = ['parcel_{}'.format(int(i)) for i in component_in]
            component_map  = {i: x for i, x in enumerate(component_in)}

            sub_bulk_file = bulk_files[0]
            for sub_bulk_file in bulk_files:
                id = sub_bulk_file.split('/')[-1].split('_')[0]

                empty_matrix = np.zeros((len(component_in), len(component_in)))
                idxs   = np.triu_indices(len(component_in), k=1)
                dat_in = np.loadtxt(sub_bulk_file)
                bulk_name = bulk_name.replace('rfmri_', 'rest_funccon_')

                edge_list = ['{}_p{}_p{}'.format(bulk_id, component_in[x], component_in[y]) for x, y in
                             zip(idxs[0], idxs[1])]
                subj_data = pd.DataFrame(dat_in).transpose()
                subj_data.columns = edge_list

                subj_data.insert(0, 'id', id)
                data_list.append(subj_data)

        out_file = os.path.join(cur_bulk_dir, '{}_dataframe.csv'.format(bulk_id))
        print('Writing:   {}'.format(out_file))
        field_df = pd.concat(data_list)
        field_df.to_csv(out_file, index=None)
        print(len(bulk_files))

    old_code = False
    if old_code == True:
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


