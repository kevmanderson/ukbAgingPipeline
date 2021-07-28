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
import textwrap
from functools import reduce
from argparse import RawTextHelpFormatter
import utilities.utilities as utilities
log = logging.getLogger('ukb')

def main():
    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    log.info('RUNNING: STAGE = RUN_PHESANT')


    parser = argparse.ArgumentParser(description=textwrap.dedent('''
                         Prepare data for phenome-wide association scan. 
                        ''')
                        )
    parser.add_argument('--phesant_csv', dest='phesant_csv', required=True, action='append', nargs='+', default=None,
                        help=textwrap.dedent('''
                         Repeatable. Full filepath to UKB csv data files produced in previous steps. 
                         This could be *enc_ukb file that has been converted to csv format.
                         Or it could be bulk data (e.g. RSFA) previously formatted for csv.
                        ''')
                        )
    parser.add_argument('--write_metadata', dest='write_metadata', required=False, default=None,
                        help=textwrap.fill('''
                            ''')
                        )
    parser.add_argument('--combined_metadata_file', dest='combined_metadata_file', required=False, default=None,
                        help = 'Path for for where to write the combined and formatted metadata file.'
                        )
    parser.add_argument('--phesant_visits', dest='phesant_visits', required=False, default='0;1;2;3',
                        help='UKB assessment visits to be included in the PHESANT analysis.'
                        )
    parser.add_argument('--prepare_data', dest='prepare_data', required=False, action='store_true', default=False,
                        help='Binary flag. If true, combine and format data for input to PHESANT pipeline'
                        )
    parser.add_argument('--run_phesant', dest='run_phesant', required=False, action='store_true', default=False,
                        help = 'Binary flag. If true, run PHESANT pipeline.'
                        )
    parser.add_argument('--outcome_info', dest='outcome_info', required=True, default=False,
                        help = 'PHESANT style "outcome-info.tsv" file with UKB metadata information.'
                        )
    parser.add_argument('--field_file', dest='field_file', required=True, default=False,
                        help = 'UKB field information. Obtained from here: https://biobank.ndph.ox.ac.uk/showcase/schema.cgi'
                        )
    parser.add_argument('--showcase', dest='showcase', required=True, default=False,
                        help = 'Additional UKB metadata information from here: https://biobank.ndph.ox.ac.uk/showcase/schema.cgi'
                        )
    parser.add_argument('--nparts', dest='nparts', required=False, default=30, help='Number of data splits for parallel PHESANT run')
    parser.add_argument('--out_dir', dest='out_dir', required=True, help="Directory to write PHESANT formatted csv files")
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    args = parser.parse_args()

    tmp = True
    if tmp == False:
        parser = argparse.ArgumentParser()
        args   = parser.parse_args()
        args.out_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw'
        args.phesant_csv = ['/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_3.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_2.csv',
                            '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_3.csv']
        args.field_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/field.txt'
        args.showcase  = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/Data_Dictionary_Showcase.tsv'
        args.outcome_info = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/outcome-info.tsv'
        args.write_metadata = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/aging-outcome-info.tsv'
        args.phesant_visits   = '0;1;2;3'
        args.phesant_visits = args.phesant_visits.split(';')
        args.nparts = '30'

    # read metadata info for PHESANT
    outcome_df    = pd.read_table(args.outcome_info)
    ukb_field_df  = pd.read_table(args.field_file)
    showcase_df   = pd.read_table(args.showcase)

    # first merge field/showcase metadata
    meta_df = ukb_field_df.merge(showcase_df, left_on='field_id', right_on='FieldID')

    # allow genetic and some related fields to run through PHESANT pipelines (don't exclude them like normal PHESANT run)
    gene_cols          = outcome_df.FieldID[outcome_df['EXCLUDED'] == 'YES-GENETIC'].tolist()
    dont_exclude_these = [int(x) for x in [26500, 31, 21003, 54] + gene_cols]

    # the EXCLUDED column in "outcome_df" acts as a data filter within PHESANT.
    # Set the exclusion criteria to nan for the fields defined above
    outcome_df.loc[outcome_df.FieldID.isin(dont_exclude_these), 'EXCLUDED'] = np.nan

    # for now, just examine variables that were collected in person
    # i.e., exclude online collections where we would have to hand calculate age
    # see here: https://biobank.ndph.ox.ac.uk/showcase/schema.cgi?id=9
    meta_df = meta_df.loc[~meta_df['Path'].str.contains('Online follow-up')]


    # We have to organize the data a bit in order to run through PHESANT properly
    # specifically, we split data into their associated UKB visits 

    # read each basket of data and split by visit
    # --------------
    visit_arr = args.phesant_visits
    step = 10000
    df_dict = {}
    for ukb_csv in args.phesant_csv:
        if type(ukb_csv) == list:
            ukb_csv = ukb_csv[0]
        log.info('Reading: {}'.format(ukb_csv))
        bucket_id, visit_dfs = read_and_reduce_csv(ukb_csv=ukb_csv,
                                                    meta_df=meta_df,
                                                    step=step,
                                                    visit_arr=visit_arr)
        df_dict[bucket_id] = visit_dfs

    # merge data from each input csv file (keep separated by visit)
    # --------------
    merge_dict = {}
    for visit in args.phesant_visits:
        print(visit)

        # get rid of any empty dataframes, then merge
        merge_dfs = [df_dict[key][visit] for key in df_dict.keys() if df_dict[key][visit].shape[0] != 0]

        # lambda magic
        df = reduce(lambda df1, df2: pd.merge(df1, df2, on='xeid', how='outer', suffixes=('', '_rm')), merge_dfs)

        # reduce duplicated columns, which can occur if the same data field was present in multiple bulk csv inputs
        keep_cols         = df.columns[~df.columns.str.contains('_rm')].tolist()

        # write combined csvs to file
        merge_dict[visit] = df[keep_cols]
        phesant_write     = os.path.join(args.out_dir, 'phesant_visit{}.csv'.format(visit))
        merge_dict[visit].to_csv(phesant_write, index=None)


    # add metadata to PHESANT if we're dealing with RSFA and RSFC data
    # --------------
    for bulk_id in ['25750', '25751', '25752', '25753', '25754', '25755']:
        if any(bulk_id in x for x in df_dict.keys()):
            meta_df         = pd.read_csv(os.path.join(args.out_dir, 'bulk_{}_metadata.csv'.format(bulk_id)))
            meta_df.columns = meta_df.columns.str.replace('field_id', 'FieldID')
            # add these columns, populated with nan vals
            nan_fields = ['TRAIT_OF_INTEREST',
                            'EXCLUDED',
                            'CAT_MULT_INDICATOR_FIELDS',
                            'CAT_SINGLE_TO_CAT_MULT',
                            'DATE_CONVERT',
                            'DATA_CODING']
            row = pd.Series({'FieldID': bulk_id,
                                'Path': meta_df['path'][0],
                                'ValueType':'Continuous',
                                'Category': meta_df['category'][0],
                                'Field': meta_df.loc[0]['field']})
            for f in nan_fields: row[f]=np.nan
            row_series = pd.Series(row)
            outcome_df = outcome_df.append(row_series, ignore_index=True)

    # write this metadata file to be used with PHESANT
    # --------------

    # get rid of codes for antibiotics taken in the last 3 months (i.e. 20199, very polymorphic and uninteresting)
    outcome_df.loc[outcome_df.FieldID == 20199, 'CAT_MULT_INDICATOR_FIELDS'] = np.nan
    outcome_df.to_csv(args.write_metadata, sep='\t', index=None)

    # write visit specific phesant formatted csv files
    # --------------
    merge_dict = {}
    for visit in args.phesant_visits:
        print(visit)
        merge_dict[visit] = pd.read_csv(os.path.join(args.out_dir, 'phesant_visit{}.csv'.format(visit)), low_memory=False)

    # PHESANT covariates
    # --------------
    covar_sex = merge_dict['0'][['xeid', 'x31_0_0']]
    for visit in args.phesant_visits:
        covar_df = covar_sex.merge(merge_dict[visit][['xeid', 'x54_{}_0'.format(visit)]], on='xeid')
        covar_df.to_csv(os.path.join(args.out_dir, 'phesant_covar_vis{}.csv'.format(visit)), index=None)


    # create a dataframe with common covariates for faster loading (particularly genetic PCs, which take a long time to load from SQL)
    pull_covars = [
        'x22001_',
        'x31_',
        'x22006_',
        'x22019_',
        'x22020_',
        'x22027_',
        'x26755_',
        'x26856_',
        'x26721_',
        'x26822_',
        'x26537_',
        'x25734_',
        'x25735_',
        'x25756_',
        'x25757_',
        'x25758_',
        'x26521_',
        'x25746_',
        'x25741_',
        'x25742_',
        'x25744_',
        'x21001_',
        'x22009_'
    ]
    common_covar_list = []
    for visit in args.phesant_visits:
        match_cols = []
        for substr in pull_covars:
            match_idxs = merge_dict[visit].columns.str.contains(substr)
            match_cols.append(merge_dict[visit].columns[match_idxs].tolist())
        flat_list = [item for sublist in match_cols for item in sublist]
        covar_df = merge_dict[visit][['xeid'] + flat_list]
        common_covar_list.append(covar_df)

    common_covar_df = reduce(lambda df1, df2: pd.merge(df1, df2, on='xeid', how='outer', suffixes=('', '_rm')), common_covar_list)
    # write csv
    common_covar_df.to_csv(os.path.join(args.out_dir, 'phesant_common_covars.csv'), index=None)
    #x = pd.read_csv(os.path.join(args.out_dir, 'phesant_common_covars.csv'))
    
    # write hd5
    common_covar_df.to_hdf(os.path.join(args.out_dir, 'phesant_common_covars.h5'), 'xeid', mode='w')
    #z = pd.read_hdf(os.path.join(args.out_dir, 'phesant_common_covars.h5'), 'xeid')


def read_and_reduce_csv(ukb_csv, meta_df, step, visit_arr):
    '''


    Parameters
    ----------
    ukb_csv: full path to csv file
    meta_df: dataframe with UKB metadata
    step:    number of rows to read when chunking through csv file
    visit_arr: array of UKB visits (0-3) to include

    Returns
    -------

    '''

    # name of the csv file
    bucket_id = ukb_csv.split('/')[-1].replace('.csv', '')

    # read csv header (w/ column names). Put them into format for PHESANT header style (i.e., x{}_{}_{})
    hdr         = pd.read_csv(ukb_csv, nrows=1)
    hdr.columns = ['x{}'.format(x) for x in hdr.columns.str.replace('-', '_').str.replace('[.]', '_')]

    # the input dataframes contain more data than will make it into the final database
    # we want to pull fields for which we have metadata, then split them by UKB assessment center visit
    # b/c sample size differs for each visit, we eliminate a lot of null values by doing these visit splits
    visitcol_dict = {}
    grep_fields   = ['x{}_'.format(x) for x in meta_df.field_id]
    for visit in visit_arr:
        visitcol_dict[visit] = {}
        this_visit_cols = hdr.columns[hdr.columns.str.contains('_{}_'.format(visit))]
        use_cols = [x for x in this_visit_cols if '{}_'.format(x.split('_')[0]) in grep_fields]
        visitcol_dict[visit] = use_cols


    # read the input csv in chunks
    # get rid of subjects who have full missingness
    visit_dfs = {}
    reader = pd.read_csv(ukb_csv, chunksize=step, low_memory=False, encoding="ISO-8859-1", )
    for i, ukb_df in enumerate(reader):
        print(i)
        log.info('Reading Chunk-size: {}, Chunk: {}'.format(step, i))

        # convert columns to PHESANT format (i.e., x{}_{}_{})
        ukb_df.columns = ['x{}'.format(x) for x in ukb_df.columns]
        ukb_df.columns = ukb_df.columns.str.replace('-', '_').str.replace('[.]', '_')

        # separate data by visits
        for visit in visit_arr:
            if i == 0:
                visit_dfs[visit] = []
            #
            columns_to_write = visitcol_dict[visit]
            vis_cols      = list(ukb_df.columns[ukb_df.columns.str.contains('_{}_'.format(visit))])

            # get rid of gene metadata columns (no actual data)
            exclude_these = ukb_df.columns[ukb_df.iloc[0] == '#ukbgene'].tolist()
            columns_to_write  = [x for x in columns_to_write if x not in exclude_these]

            # pull columns for the given visit
            vis_df     = ukb_df[['xeid'] + vis_cols]
            check_cols = [x for x in columns_to_write if x != 'xeid']

            # get rid of subs with all NAs (e.g. subjects that never had MRI/came in for visit 2)
            keep_rows  = vis_df[check_cols].isna().sum(1) != len(check_cols)
            write_df   = vis_df[['xeid'] + columns_to_write]
            write_df   = write_df.loc[keep_rows]
            visit_dfs[visit].append(write_df)

    # after reading, concatenate all of the chunks into single dataframes (split by visit)
    for visit in visit_arr:
        if type(visit_dfs[visit]) != pd.DataFrame: 
            visit_dfs[visit] = pd.concat(visit_dfs[visit])

    return bucket_id, visit_dfs




if __name__ == "__main__":
    main()
