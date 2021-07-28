#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# Convert decrypted UKB data to readable formats

import os
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
    This is a wrapper for the "ukbconv" utility
    Detailed info here: https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide

    This is called within in a python script to make sure the data end up in the appropriate directories.

    Parameters
    ----------
    enc_ukb: string, required
        full path to the decrypted UKB .enc_ukb file
        e.g., /Users/me/ukb12345.enc_ukb

    Returns
    -------

    '''
    logging.basicConfig(level=logging.DEBUG)
    utilities.configLogging()
    log.info('buckner-lab-ukb-pipeline')
    log.info('RUNNING: CONVERT UKB')

    parser = argparse.ArgumentParser()
    parser.add_argument('--enc_ukb', dest='enc_ukb', required=True, action='append', nargs='+', default=None)
    parser.add_argument('--ukb_conv', dest='ukb_conv', required=False, default=None)
    parser.add_argument('--fields_to_convert', dest='fields_to_convert', required=False, default=None)
    parser.add_argument('--bulk_field', dest='bulk_field', required=False, action='append', nargs='+', default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True)
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    args = parser.parse_args()
    print(args)

    tmp = True
    if tmp == False:
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.enc_ukb = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc_ukb']
        #args.category_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/funpack/categories.tsv'
        args.out_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging'
        args.fields_to_convert = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/fields_to_convert.txt'
        args.showcase_df = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/Data_Dictionary_Showcase.tsv'
    


    # ukb utility to convert encrypted data
    #ukbconv   = os.path.join(args.out_dir, 'data/ukb/external/ukbconv')
    #raw_dir   = os.path.join(args.out_dir, 'data/ukb/raw')
    write_dir = os.path.join(args.out_dir, 'slurm')

    # read UKB fields to extract
    # this is formatted as described here: https://git.fmrib.ox.ac.uk/fsl/funpack
    tmp = False
    if tmp == None:

        # write list of fields to convert
        #field_arr  = utilities.read_funpack_categories(args.category_file)
        #field_arr_series = pd.Series(field_arr)
        #nomatch = field_arr_series[~pd.Series(field_arr).isin(showcase_df['FieldID'])]

        showcase_df  = pd.read_table('/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/Data_Dictionary_Showcase.tsv')
        
        # PHESANT has lots of useful information about fields to exclude
        phesant_df   = pd.read_table(args.outcome_info)
        # exclude phenos within these cats
        exclude_cats = [
            'YES-CAT-SIN-MUL-VAL',
            'YES-DATE',
            'YES-LOCATION',
            'YES-POLYMORPHIC',
            'YES-PROCESSING',
            'YES-SENSITIVE',
            'YES-SIMONS',
            'YES-SUPERSEDED'
        ]
        phesant_filt_df = phesant_df.loc[~phesant_df['EXCLUDED'].isin(exclude_cats)]
        showcase_df = showcase_df.loc[showcase_df['FieldID'].isin(phesant_filt_df['FieldID'])]
        download_df = showcase_df[['FieldID', 'Field']]
        download_df.to_csv('/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/ukb_download_df.txt', sep='\t', index=None)

        field_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/fields_to_convert.txt'
        field_arr  = list(np.sort(list(set(download_df['FieldID']))))
        with open(field_file, 'w') as f:
            [f.writelines('{}\n'.format(field) for field in field_arr)]
    else:
        # Convert all fields
        field_file = None

    # decrypt all input enc/key pairs
    for ukb_enc in args.enc_ukb:
        for output_format in ['csv', 'r']:
            o_name      = ukb_enc[0].split('/')[-1].replace('.enc_ukb', '')
            base_cmd = f"cd {args.out_dir}\n"
            if field_file is not None:
                convert_cmd = f"{base_cmd}{args.ukb_conv} {ukb_enc[0]} {output_format} -o./{o_name} -i./{field_file}"
            else:
                convert_cmd = f"{base_cmd}{args.ukb_conv} {ukb_enc[0]} {output_format} -o./{o_name}"

            write_file = os.path.join(args.out_dir, 'convert_{}_{}.bash'.format(o_name, output_format))
            utilities.write_cmd(cmd=convert_cmd, write_file=write_file)



if __name__ == "__main__":
    main()



