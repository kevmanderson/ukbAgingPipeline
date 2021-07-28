#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# Convert decrypted UKB data to readable formats

import os
import glob
import datetime
import calendar
import argparse
import logging
import numpy as np
import pandas as pd
import utilities.utilities as utilities

log = logging.getLogger('ukb')

def main():
    '''
    Concat all bulk data into a single csv

    Parameters
    ----------


    Returns
    -------

    '''

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()

    parser = argparse.ArgumentParser()
    parser.add_argument('--ukb_fetch', dest='ukb_fetch', required=True, default=None)
    parser.add_argument('--ref_dir', dest='ref_dir', required=True, default=None)
    parser.add_argument('--bulk_field', dest='bulk_field', required=True, action='append', nargs='+', default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True)
    args = parser.parse_args()
    print(args)

    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    args.bulk_field = [['rfmri_full_25:25750', 'rfmri_full_100:25751']]
    args.out_dir    = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw'
    args.ukb_fetch  = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/ukb_tools/ukbfetch'
    args.slurm_dir  = '/gpfs/milgram/project/holmes/kma52/buckner_aging/slurm'
    args.enc_ukb    = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc_ukb'
    args.ref_dir    = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files'
    log.info('RUNNING: STAGE = BULK_TO_CSV')

    # point to directories
    bulk_dir  = os.path.join(config_json['base_dir'], 'data/ukb/bulk')
    raw_dir   = os.path.join(config_json['base_dir'], 'data/ukb/raw')

    # iterate through all of the bulk fields to download
    for field in args.bulk_field:
        bulk_name  = field[0].split(':')[0]
        bulk_id    = field[0].split(':')[1]
        bulk_dir   = os.path.join(args.out_dir, bulk_name)
        for visit in ['2','3']:
            read_restbulk_and_write(bulk_dir, args.ref_dir, args.out_dir, bulk_id, visit)


def read_restbulk_and_write(bulk_dir, ref_dir, out_dir, bulk_id, visit):

    # get the descriptive bulk_name and ukb bulk_id
    bulk_files = glob.glob(os.path.join(bulk_dir, '*{}*_{}_0.txt'.format(bulk_id, visit)))

    # 25 component features
    # read the mapping between index and ROI number
    if bulk_id in ['25750', '25752', '25754']:
        n_comp     = 21
        total_comp = 25
        comp_path  = os.path.join(ref_dir, 'rfMRI_GoodComponents_d25_v1.txt')
        log.debug('Read component feature names from FILE ("{}")'.format(comp_path))
        with open(comp_path, 'r') as f: comp_names = f.readline().split()

    # 100 component features
    # read the mapping between index and ROI number
    elif bulk_id in ['25751', '25753', '25755']:
        n_comp     = 55
        total_comp = 100
        comp_path  = os.path.join(ref_dir, 'rfMRI_GoodComponents_d100_v1.txt')
        log.debug('Read component feature names from FILE ("{}")'.format(comp_path))
        with open(comp_path, 'r') as f: comp_names = f.readline().split()

    # RSFA
    if bulk_id in ['25754', '25755']:
        log.debug('RSFA')
        n_values  = n_comp
        col_names = ['{}-{}.{}'.format(bulk_id, visit, i) for i in range(len(comp_names))]
        long_col_names = ['MRI Rest RSFA region-{}, {} components'.format(x, total_comp) for x in comp_names]

    #'111'
    # func-con
    elif bulk_id in ['25750', '25751', '25752', '25753']:
        log.debug('RSFC')
        n_values = ((n_comp * n_comp) - n_comp) / 2
        zero_df  = pd.DataFrame(np.zeros((n_comp, n_comp)))
        zero_df.columns = comp_names
        zero_df.index   = comp_names
        zero_df = zero_df.where(np.triu(np.ones(zero_df.shape), k=1).astype(np.bool))
        zero_df = zero_df.stack().reset_index()
        edges   = zip(zero_df[zero_df.columns[0]], zero_df[zero_df.columns[1]])
        col_names = ['{}-{}.{}'.format(bulk_id, visit, x) for x in range(len(zero_df[zero_df.columns[0]]))]
        long_col_names = ['MRI Rest RSFC region-{} to region-{}, {} components'.format(x[0], x[1], total_comp) for x in edges]


    # make/write metadata
    # ---------
    meta_df = pd.DataFrame({'field_id': bulk_id,
                            'field': long_col_names,
                            'col_name': col_names,
                            'category': 111,
                            'visit': visit,
                            'path': 'UK Biobank Assessment Centre > Imaging > Brain MRI > Resting functional brain MRI',
                            'valuetype': 'Continuous'})
    csv_path = os.path.join(out_dir, 'bulk_{}_metadata.csv'.format(bulk_id, visit))
    meta_df.to_csv(csv_path, index=None)


    # read the individual data files
    # ---------
    data_list  = []
    for i,bulk_path in enumerate(bulk_files):
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
    bulk_df.columns = ['eid'] + col_names

    # write bulk data to csv
    csv_path = os.path.join(raw_dir, 'bulk_{}_{}.csv'.format(bulk_id, visit))
    log.info('Writing bulk data to FILE ("{}")'.format(csv_path))
    bulk_df.to_csv(csv_path, index=None)



if __name__ == "__main__":
    main()





