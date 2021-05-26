
import os
import sys
import wget
import glob
import numpy as np
import pandas as pd
import logging

from utilities.utilities import make_symlink
log = logging.getLogger('ukb')


def read_and_reduce_csv(ukb_csv, meta_df, step):
    bucket_id = ukb_csv.split('/')[-1].replace('.csv', '')

    # read csv header
    hdr         = pd.read_csv(ukb_csv, nrows=1)
    hdr.columns = ['x{}'.format(x) for x in hdr.columns.str.replace('-', '_').str.replace('[.]', '_')]

    # identify columns to pull from the larger dataframe
    visitcol_dict = {}
    grep_fields   = ['x{}_'.format(x) for x in meta_df.field_id]
    for visit in ['0','1','2','3']:
        visitcol_dict[visit] = {}
        this_visit_cols = hdr.columns[hdr.columns.str.contains('_{}_'.format(visit))]
        use_cols = [x for x in this_visit_cols if '{}_'.format(x.split('_')[0]) in grep_fields]
        visitcol_dict[visit] = use_cols

    visit_dfs = {}
    reader = pd.read_csv(ukb_csv, chunksize=step, low_memory=False, encoding="ISO-8859-1", )
    for i, ukb_df in enumerate(reader):
        log.debug('Reading Chunk-size: {}, Chunk: {}'.format(step, i))
        ukb_df.columns = ['x{}'.format(x) for x in ukb_df.columns]
        ukb_df.columns = ukb_df.columns.str.replace('-', '_').str.replace('[.]', '_')

        for visit in ['0','1','2','3']:
            if i == 0:
                visit_dfs[visit] = []
            columns_to_write = visitcol_dict[visit]
            vis_cols   = list(ukb_df.columns[ukb_df.columns.str.contains('_{}_'.format(visit))])
            exclude_these = ukb_df.columns[ukb_df.iloc[0] == '#ukbgene'].tolist()
            columns_to_write   = [x for x in columns_to_write if x not in exclude_these]
            vis_df     = ukb_df[['xeid'] + vis_cols]
            check_cols = [x for x in columns_to_write if x != 'xeid']
            keep_rows  = vis_df[check_cols].isna().sum(1) != len(check_cols)
            write_df   = vis_df[['xeid'] + columns_to_write]
            write_df   = write_df.loc[keep_rows]
            visit_dfs[visit].append(write_df)

    for visit in ['0','1','2','3']:
        visit_dfs[visit] = pd.concat(visit_dfs[visit])

    return bucket_id, visit_dfs


def prep_data_for_phesant(config_json, args):
    '''Concat all bulk data into a single csv'''

    log.info('RUNNING: STAGE = PREP_DATA_FOR_PHESANT')

    # ukb showcase metadata
    ukb_field_df = pd.read_table(os.path.join(config_json['repo_dir'], 'ref_files/field.txt'))
    showcase_df  = pd.read_csv(os.path.join(config_json['repo_dir'], 'ref_files/showcase.csv'))
    meta_df      = ukb_field_df.merge(showcase_df, left_on='field_id', right_on='FieldID')
    outcome_df   = pd.read_table(os.path.join(config_json['repo_dir'], 'ref_files/outcome-info.tsv'))

    # allow genetic and some related fields to run through PHESANT pipelines (don't exclude)
    gene_cols = [str(x) for x in outcome_df.FieldID[outcome_df['EXCLUDED'] == 'YES-GENETIC'].tolist()]
    dont_exclude_these = [int(x) for x in ['26500', '31', '21003', '54'] + gene_cols]
    outcome_df.loc[outcome_df.FieldID.isin(dont_exclude_these), 'EXCLUDED'] = np.nan

    # for now, just examine variables that were collected in person
    # exclude online collections where we have to hand calculate age
    # see here: https://biobank.ndph.ox.ac.uk/showcase/schema.cgi?id=9
    meta_df = meta_df.loc[~meta_df['Path'].str.contains('Online follow-up')]

    # list of csvs to process/combined
    ukb_csv_list = [x['ukb_enc'].replace('.enc', '.csv') for i, x in enumerate(config_json['ukb_encs'])]
    if type(args.phesant_csv) == str: args.phesant_csv=[args.phesant_csv]
    ukb_csv_list = ukb_csv_list + args.phesant_csv

    # read each basket of data and split by visit
    step = 10000
    df_dict = {}
    for ukb_csv in ukb_csv_list:
        bucket_id, visit_dfs = read_and_reduce_csv(ukb_csv, meta_df, step)
        df_dict[bucket_id] = visit_dfs

    # merge by visit
    merge_dict = {}
    for visit in ['0','1','2','3']:
        print(visit)
        # get rid of any empty dataframes, then merge
        merge_dfs = [df_dict[key][visit] for key in df_dict.keys() if df_dict[key][visit].shape[0] != 0]
        df        = reduce(lambda df1, df2: pd.merge(df1, df2, on='xeid', how='outer', suffixes=('', '_rm')), merge_dfs)
        keep_cols = df.columns[~df.columns.str.contains('_rm')].tolist()
        merge_dict[visit] = df[keep_cols]
        phesant_write = os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_visit{}.csv'.format(visit))
        merge_dict[visit].to_csv(phesant_write, index=None)

    # add metadata to PHESANT if we're dealing with RSFA and RSFC data
    for bulk_id in ['25750', '25751', '25752', '25753', '25754', '25755']:
        if any(bulk_id in x for x in df_dict.keys()):
            meta_df = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/bulk_{}_metadata.csv'.format(bulk_id)))
            meta_df.columns = meta_df.columns.str.replace('field_id', 'FieldID')
            meta_df.columns = meta_df.columns.str.replace('field_id', 'FieldID')
            nan_fields = ['TRAIT_OF_INTEREST', 'EXCLUDED', 'CAT_MULT_INDICATOR_FIELDS', 'CAT_SINGLE_TO_CAT_MULT', 'DATE_CONVERT', 'DATA_CODING', ]
            row = pd.Series({'FieldID': bulk_id, 'Path': meta_df['path'][0], 'ValueType':'Continuous',
                             'Category': meta_df['category'][0], 'Field': meta_df.loc[0]['field']})
            for f in nan_fields: row[f]=np.nan
            row_series = pd.Series(row)
            outcome_df = outcome_df.append(row_series, ignore_index=True)

    # write this metadata file to be used with PHESANT
    outcome_df.loc[outcome_df.FieldID == 20199, 'CAT_MULT_INDICATOR_FIELDS'] = np.nan
    outcome_df.to_csv(os.path.join(config_json['repo_dir'], 'ref_files/aging-outcome-info.tsv'), sep='\t', index=None)

    # create covariate dataframe
    merge_dict = {}
    for visit in ['0','1','2','3']:
        print(visit)
        merge_dict[visit] = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_visit{}.csv'.format(visit)), low_memory=False)

    # PHESANT covariates
    covar_sex = merge_dict['0'][['xeid', 'x31_0_0']]
    for visit in ['0', '1', '2', '3']:
        covar_df = covar_sex.merge(merge_dict[visit][['xeid', 'x54_{}_0'.format(visit)]], on='xeid')
        covar_df.to_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_covar_vis{}.csv'.format(visit)), index=None)





def run_phesant(config_json, args):
    log.info('RUNNING: STAGE = RUN_PHESANT')
    log.debug(args)

    write_dir = os.path.join(config_json['base_dir'], 'slurm')

    # age_var     = args.age_var #'x21003_2_0'
    visit_arr = args.phesant_visits.split(';')

    # set up filepaths
    phesant_dir      = os.path.join(config_json['repo_dir'], 'external/PHESANT/WAS')
    resDir           = os.path.join(config_json['base_dir'], 'data/ukb/phesant/')
    variablelistfile = os.path.join(config_json['repo_dir'], 'ref_files/aging-outcome-info.tsv')
    datacodingfile   = os.path.join(config_json['repo_dir'], 'ref_files/data-coding-ordinal-info.txt')

    nparts = 30
    for visit in visit_arr:
        age_var        = 'x21003_{}_0'.format(visit)
        confounderfile = os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_covar_vis{}.csv'.format(visit))
        phenofile      = os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_visit{}.csv'.format(visit))

        for part in range(1,nparts+1):
            phesant_cmd = f'''
cd {phesant_dir}
Rscript phenomeScan.r \\
        --phenofile={phenofile} \\
        --variablelistfile={variablelistfile} \\
        --confounderfile={confounderfile} \\
        --datacodingfile={datacodingfile} \\
        --traitofinterest="{age_var}" \\
        --resDir={resDir} \\
        --userId="xeid" \\
        --sensitivity \\
        --numParts={nparts} --partIdx={part} \\
        --ageistraitofinterest \\
        --genetic="FALSE" \\
        --standardise="TRUE" \\
        --visit={visit} \\
        --file_prepend=phewas_visit{visit} \\
        --save'''

            # write executable to file
            if args.slurm == True:
                log.info('submit slurm')
                write_file = os.path.join(write_dir, 'phesant_visit{}_part{}-{}'.format(visit, part, nparts))
                utilities.submit_slurm(
                    utilities.write_slurm(slurm_file=write_file,
                                          partition=args.slurm_partition, cmd=phesant_cmd,
                                          jobName='{}_{}'.format(visit, part), stime='6:00:00', nthreads=12))
            else:
                write_file = os.path.join(write_dir, 'phesant_visit{}_part{}-{}'.format(visit, part, nparts))
                utilities.write_cmd(cmd=phesant_cmd, write_file=write_file)


















