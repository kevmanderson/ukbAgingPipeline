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
    Before downloading bulk, data first need to make lists of subjecst that possess data

    Parameters
    ----------


    Returns
    -------

    '''

    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    log.info('RUNNING: STAGE = RUN_PHESANT')


    # set up filepaths/directories
    sql_dir       = os.path.join(config_json['base_dir'], 'data/ukb/sql')

    # read variable information
    ukb_meta      = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv'))
    ordinal_codes = pd.read_csv(os.path.join(config_json['repo_dir'], 'ref_files/data-coding-ordinal-info.txt'))
    codings       = pd.read_csv(os.path.join(config_json['repo_dir'], 'ref_files/codings.csv'))

    # combined PHESANT regression results into single table
    phesant_df = combine_phewas_results(config_json, args)
    phesant_df.to_csv(os.path.join(config_json['base_dir'], 'output/sqlite_files/phesant_stats.csv'), index=None)

    # define fields to not include in db ( plus icd9+icd10 codes)
    exclude_cats   = ['YES-ACE', 'YES-CAT-SIN-MUL-VAL', 'YES-DATE', 'YES-POLYMORPHIC', 'YES-SENSITIVE']
    exclude_fields = ukb_meta.loc[ukb_meta['excluded'].isin(exclude_cats), 'fieldid'].tolist() + [41202, 41203, 41204, 41205, 41270, 41271, 41216]

    # make a list of all column names to commit to sql
    ukb_csv_list = []
    for visit in range(0,4):
        ukb_csv_list.append(os.path.join(config_json['base_dir'], 'data/ukb/raw', 'phesant_visit{}.csv'.format(visit)))

    # read headers from each datafile
    hdr_cols = read_csv_headers(ukb_csv_list)

    # bin each column as a str, int, real, or date datatype
    table_types = dict(zip(['str','int','real','datetime'], ['VARCHAR', 'INTEGER', 'REAL', 'REAL']))
    table_cols  = field_to_tables(table_types=table_types, ukb_meta=ukb_meta, ukb_cols=hdr_cols, exclude_fields=exclude_fields)

    # --------
    # METADATA
    # --------
    # create a metadata dataframe of all fields in the SQL db
    sql_meta_df = make_sql_df(table_cols)

    # subset ukb metadata to include only data going into the SQL db
    ukb_meta    = ukb_meta.loc[ukb_meta.fieldid.isin(sql_meta_df.fieldid)]

    # create metadata for categorical multiple fields
    multi_metadata_df = make_cat_multi_metadata(ukb_meta, sql_meta_df)

    # merge everything into a final metadata df to commit to the database
    write_meta_df = merge_metadata(sql_meta_df, multi_metadata_df, ukb_meta)
    #ukb_meta      = pd.read_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv'))

    # --------
    # DATA
    # --------

    import sqlite3
    sql_path = os.path.join(sql_dir, 'ukb_db.sqlite')
    sql_path = os.path.join(sql_dir, 'ukb_db_smaller.sqlite')
    con = sqlite3.connect(database=sql_path)

    # initialize the empty tables in the sqlite db
    sql_resp = con.executescript("".join(initiate_tables(table_types)))

    # commit the large UKB to SQL file in chunks
    step     = 10000
    for ukb_csv in ukb_csv_list:
        reader   = pd.read_csv(ukb_csv, chunksize=step, low_memory=False, encoding="ISO-8859-1",)
        for i, ukb_chunk in enumerate(reader):
            log.debug('Reading Chunk-size: {}, Chunk: {}'.format(step, i))
            chunk_cols = set(ukb_chunk.columns) - set('xeid')
            for table_name, sql_type in table_types.items():
                table_fields = table_cols[table_name]
                sql_insert   = con.executemany(*insert_data_to_sql(ukb_chunk, table_name, table_fields))

    #
    con.commit()

    # con.execute('DROP INDEX str_index')
    # con.execute('DROP INDEX int_index')
    # con.execute('DROP INDEX real_index')
    # con.execute('DROP INDEX dt_index')

    # create index for each table (increases size a ton, but speeds up searches)
    con.execute('CREATE INDEX str_index ON str (field, time, array)')
    con.execute('CREATE INDEX int_index ON int (field, time, array)')
    con.execute('CREATE INDEX real_index ON real (field, time, array)')
    con.execute('CREATE INDEX dt_index ON datetime (field, time, array)')


    write_meta_df = add_covariate_presets(write_meta_df)
    write_meta_df['pheno_category'] = [x.split(' > ')[-1] for x in write_meta_df['path']]
    fieldid_arr = write_meta_df.loc[write_meta_df['col_name'].isnull(), 'fieldid'].astype(str)
    array_arr   = write_meta_df.loc[write_meta_df['col_name'].isnull(), 'array'].astype(str)
    visit_arr   = write_meta_df.loc[write_meta_df['col_name'].isnull(), 'visit'].astype(int).astype(str)

    fill_col_name = 'x' + fieldid_arr + '_' + visit_arr + '_' + array_arr
    write_meta_df.loc[write_meta_df['col_name'].isnull(), 'col_name'] = fill_col_name

    write_meta_df['sexed'] = write_meta_df['sexed_x']
    write_meta_df['units'] = write_meta_df['units_x']
    write_meta_df['notes'] = write_meta_df['notes_x']
    write_meta_df = write_meta_df[write_meta_df.columns[~write_meta_df.columns.str.contains('_x')]]
    write_meta_df = write_meta_df[write_meta_df.columns[~write_meta_df.columns.str.contains('_y')]]

    write_meta_df.to_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/formatted_ukb_metadata.csv'), index=None)

    write_meta_df.to_sql("metadata", con, if_exists='replace', index=False)
    codings.to_sql('datacodes', con, if_exists='replace', index=False)
    ordinal_codes.to_sql('code_ordinal', con, if_exists='replace', index=False)
    phesant_df.to_sql('phesant_stats', con, if_exists='replace', index=False)

    covar_table = ukb_covariate_table(config_json, args)
    covar_table.to_csv(os.path.join(config_json['base_dir'], 'data/ukb/raw/covariate_table.csv'), index=None)

    covar_table.to_sql('covariate_table', con, if_exists='replace', index=False)

    con.commit()


def add_covariate_presets(write_meta_df):
    write_meta_df['covar_type'] = 'centerSex'

    # diffusions/task/t2/swi fields all control for t1-eTIV
    write_meta_df.loc[write_meta_df.path.str.contains('Diffusion brain MRI'), 'covar_type'] = 'centerSexHeadXYZ_t1eTIV_dMRI'
    write_meta_df.loc[write_meta_df.path.str.contains('Task functional brain MRI'), 'covar_type'] = 'centerSexHeadXYZ_t1eTIV_task'
    write_meta_df.loc[write_meta_df.path.str.contains('T2-weighted brain MRI'), 'covar_type'] = 'centerSexHeadXYZ_t1eTIV'
    write_meta_df.loc[write_meta_df.path.str.contains('Susceptibility weighted brain MRI'), 'covar_type'] = 'centerSexHeadXYZ_t1eTIV'

    # t1/freesurfer structurals require a bit more fine-grained matching to covariates
    t1_matches = write_meta_df['path'].str.contains('T1 structural brain MRI')

    vol_matches    = write_meta_df['field'].str.contains('Volume of ')
    t1_vol_matches = [all(x) for x in zip(t1_matches, vol_matches)]
    write_meta_df.loc[t1_vol_matches, 'covar_type'] = 'centerSexHeadXYZ_t1eTIV'


    # thickness/area/grey-white all don't control for eTIV
    nonVol_matches    = write_meta_df['field'].str.contains('Mean thickness of |Area of |Grey-white contrast |Mean intensity of ')
    t1_nonVol_matches = [all(x) for x in zip(t1_matches, nonVol_matches)]
    write_meta_df.loc[t1_nonVol_matches, 'covar_type'] = 'centerSexHeadXYZ_t1NoETIV'

    return write_meta_df



def ukb_covariate_table(config_json, args):

    covar_name_map = {
        'sex': '22001',
        'center': '54',
        'headPosX': '25756',
        'headPosY': '25757',
        'headPosZ': '25758',
        't1CNR': '25756',
        't1SNR': '25734',
        'eTIV': '26521',
        'rfmriMotion': '25741',
        'rfmriInvSNR': '25744',
        'dmriOutliers': '25746',
        'tfmriMotion': '25742'
    }

    covar_type_map = {
        'centerSex':['sex','center'],
        'centerSexHeadXYZ': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ'],
        'centerSexHeadXYZ_t1NoETIV': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ', 't1CNR', 't1SNR'],
        'centerSexHeadXYZ_t1eTIV': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ', 't1CNR', 't1SNR', 'eTIV'],
        'centerSexHeadXYZ_t1eTIV_rfMRI': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ', 't1CNR', 't1SNR', 'eTIV',
                                          'rfmriMotion', 'rfmriInvSNR'],
        'centerSexHeadXYZ_t1eTIV_dMRI': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ', 't1CNR', 't1SNR', 'eTIV',
                                          'dmriOutliers'],
        'centerSexHeadXYZ_t1eTIV_task': ['sex', 'center', 'headPosX', 'headPosY', 'headPosZ', 't1CNR', 't1SNR', 'eTIV',
                                         'tfmriMotion'],
    }

    # create covariate dataframe
    df_list = []
    for key in covar_type_map.keys():
        print(key)
        ukbField = [covar_name_map[var] for var in covar_type_map[key]]
        df_list.append(pd.DataFrame({'ukbField':ukbField, 'ukbVarName':covar_type_map[key], 'covarType': key}))
    covariate_df = pd.concat(df_list)

    return covariate_df



def make_metadata(config_json, args):

    repo_dir = config_json['repo_dir']

    showcase      = pd.read_csv(os.path.join(repo_dir, 'ref_files/showcase.csv'))
    codings       = pd.read_csv(os.path.join(repo_dir, 'ref_files/codings.csv'))
    outcome_info  = pd.read_table(os.path.join(repo_dir, 'ref_files/aging-outcome-info.tsv'))
    ordinal_codes = pd.read_csv(os.path.join(repo_dir, 'ref_files/data-coding-ordinal-info.txt'))
    #regr_filter   = pd.read_csv(os.path.join(repo_dir, 'ref_files/ukb_allow_regression.csv'))

    fields        = pd.read_table(os.path.join(repo_dir, 'ref_files/field.txt'))
    fields.columns = fields.columns.str.replace('field_id','fieldid')

    # outcome_info comes from PHESANT. use this as the base reference dataframe, and add from there
    outcome_info.columns = outcome_info.columns.str.lower()
    outcome_cols = [
        'fieldid',
        'excluded',
        'cat_mult_indicator_fields',
        'cat_single_to_cat_mult',
        'date_convert',
        'valuetype'
    ]

    # showcase data from ukb website; https://biobank.ndph.ox.ac.uk/showcase/schema.cgi?id=1
    showcase.columns = showcase.columns.str.lower()
    showcase_cols = ['fieldid',
                     'path',
                     'category',
                     'field',
                     'participants',
                     'items',
                     'units',
                     'sexed',
                     'coding',
                     'notes',
                     'link',
                     ]

    showcase_merge = showcase[showcase_cols]
    outcome_info_merge = outcome_info[outcome_cols]

    # merge all the data
    meta_df = fields.merge(showcase_merge, on='fieldid')
    meta_df = meta_df.merge(outcome_info_merge, on='fieldid')
    #meta_df = meta_df.merge(regr_filter[['fieldid','allow_regression']], on='fieldid')

    write_file = os.path.join(config_json['base_dir'], 'data/ukb/raw/ukb_metadata.csv')
    meta_df.to_csv(write_file, index=None)
    return meta_df


def merge_metadata(sql_meta_df, multi_metadata_df, ukb_meta):
    # merge into single dataframe, where each individual field/visit/array combination is a row
    tmp_sql_meta_df   = sql_meta_df.loc[~sql_meta_df.fieldid.isin(multi_metadata_df.fieldid)]
    ukb_base_meta     = ukb_meta.loc[~ukb_meta.fieldid.isin(multi_metadata_df.fieldid)]
    ukb_meta_long     = pd.concat([ukb_base_meta, multi_metadata_df])

    write_meta_df = tmp_sql_meta_df[['col_name','fieldid']].merge(ukb_meta_long, on='fieldid', how='outer')
    vis_fill = [x.split('_')[1] for x in write_meta_df.col_name[write_meta_df.visit.isnull()]]
    write_meta_df.visit[write_meta_df.visit.isnull()] = vis_fill
    arr_fill = [x.split('_')[2] for x in write_meta_df.col_name[write_meta_df.array.isnull()]]
    write_meta_df.array[write_meta_df.array.isnull()] = arr_fill
    return write_meta_df



def field_to_tables(table_types, ukb_meta, ukb_cols, exclude_fields):
    '''

    :cvar
    '''

    # process the phesant column names and get the ukb field ids
    #field_ids = [x.replace('x', '').split('_')[0] if 'x' in x.split('_')[0] else x for x in ukb_cols]
    #field_ids = [x.split('_')[0] for x in field_ids]
    #field_ids = [int(x.split('#')[0]) for x in field_ids]

    ukb_cols_short = ['{}_'.format(x.split('_')[0]) for x in ukb_cols]
    ukb_col_map    = dict(zip(ukb_cols, ukb_cols_short))
    ukb_cols_short = [x for x in ukb_cols_short if x != 'xeid_']

    # keys are field_ids
    # values are phesant columns
    #fieldname_dict = dict(zip(ukb_cols, field_ids))

    table_name = 'str'
    sql_type = 'VARCHAR'

    table_field_map = {}
    table_map = dict(zip(['str', 'int', 'real', 'datetime'], [['Categorical multiple','Categorical single'],
                                                              ['Integer'],
                                                              ['Continuous'],
                                                              ['Date']]))
    for table_name, sql_type in table_types.items():
        print('{}/{}'.format(table_name, sql_type))
        # fields for this datatype (from metadata)
        table_fields     = ukb_meta['fieldid'][ukb_meta['valuetype'].isin(table_map[table_name])].tolist()
        #columns_for_this_datatype = [fieldname_dict[field] for field in table_fields if field in fieldname_dict.keys()]

        # exclude variables
        table_fields = list(set(table_fields) - set(exclude_fields))

        table_fields_fmt = ['x{}_'.format(x) for x in table_fields]

        res = [list(filter(lambda x: var in x, ukb_cols)) for var in table_fields_fmt]
        save_fields = [item for sublist in res for item in sublist]

        #table_fields_fmt = list(set(table_fields_fmt).intersection(set(ukb_cols_short)))
        #save_fields = [ukb_col_map[x] for x in table_fields_fmt]
        table_field_map[table_name] = save_fields

    return table_field_map


def insert_data_to_sql(ukb_chunk, table_name, table_fields):

    #table_cols      = ['x{}_'.format(x) for x in table_fields]
    curr_tab_fields = set(ukb_chunk.columns.to_list()).intersection(set(table_fields))

    # melt into a super long table
    ukb_melt = ukb_chunk[['xeid'] + list(curr_tab_fields)].melt(id_vars='xeid', value_vars=curr_tab_fields)
    ukb_melt = ukb_melt[ukb_melt['value'].notnull()]

    # extract ukb field from column names
    ukb_melt['field'] = [x.split('_')[0].replace('x', '') for x in ukb_melt['variable']]
    ukb_melt['time']  = [x.split('_')[1] for x in ukb_melt['variable']]
    ukb_melt['array'] = [x.split('_')[2] for x in ukb_melt['variable']]
    ukb_melt.columns = ukb_melt.columns.str.replace('xeid', 'eid')

    ukb_melt = ukb_melt[["eid", "variable", "field", "time", "array", "value"]]
    ukb_melt[ukb_melt.columns[ukb_melt.columns != 'variable']]

    # return SQL insert statement and data to insert
    return (f'INSERT INTO {table_name} values({",".join("?" * len(ukb_melt.columns))})',
            ukb_melt.values.tolist())



def initiate_tables(table_types):
    sql_make_cmds = []
    for tab_name, field_type in table_types.items():
        sql_make_cmds.append(f'DROP TABLE IF EXISTS {tab_name};')
        # names of the columns that will be populated
        field_cols      = ["eid", "variable", "field", "time", "array", "value"]
        field_col_types = ["INTEGER", "VARCHAR", "VARCHAR", "INTEGER", "INTEGER", field_type]

        cols = ','.join(map(' '.join, zip(field_cols, field_col_types)))
        sql_cmd  = f"CREATE TABLE {tab_name} ({cols}) ;"
        sql_make_cmds.append(sql_cmd)
    return sql_make_cmds


def recode_values(df, reassign):
    #  multi_code_values = codings.loc[codings['Coding'] == multi_code]
    # df = multi_code_values
    # reassign = ordinal_info.reassignments
    reassign_string = list(reassign)[0]
    recode_dict     = {x.split('=')[0]: x.split('=')[1] for x in reassign_string.split('|')}
    df['Value'] = df['Value'].replace(recode_dict)
    # remove NA
    df = df.loc[df.Value != 'NA']
    df['Value'] = pd.to_numeric(df['Value'])
    df = df.loc[df.Value >= 0]
    return df

def read_csv_headers(csv_list):
    # read headers from each datafiles
    hdr_cols = []
    for ukb_csv in ukb_csv_list:
        hdr = pd.read_csv(ukb_csv, nrows=1, encoding="ISO-8859-1")
        hdr_cols = hdr_cols + hdr.columns.tolist()
    hdr_cols = list(set(hdr_cols))
    return hdr_cols

def make_sql_df(table_cols):
    all_cols = table_cols['str'] + table_cols['int'] + table_cols['real'] + table_cols['datetime']
    all_cols = list(np.sort(all_cols))

    sql_meta_df = pd.DataFrame({'col_name':all_cols})
    sql_meta_df['fieldid'] = [pd.to_numeric(x.replace('x', '').split('_')[0]) for x in sql_meta_df['col_name']]
    sql_meta_df['visit'] = [x.split('_')[1] for x in sql_meta_df['col_name']]
    return sql_meta_df

def make_cat_multi_metadata(ukb_meta, sql_meta_df):

    # recode multiple cat fields
    multi_meta_df = ukb_meta.loc[ukb_meta['valuetype'] == 'Categorical multiple']
    multi_fields  = list(set(multi_meta_df['fieldid']))
    list_of_multi_metadata = []

    # for each unique multiple categorical data field
    for multi_field in multi_fields:
        # pull metadata info
        multi_df          = multi_meta_df.loc[multi_meta_df['fieldid'] == multi_field]

        # visits in the SQL db for this field
        multi_visits      = list(set(sql_meta_df[sql_meta_df.fieldid == multi_field].visit))

        # expand out each level of datacoding
        multi_code        = list(multi_df.coding)[0]
        multi_code_values = codings.loc[codings['Coding'] == multi_code]
        multi_code_values = multi_code_values.loc[~multi_code_values.Value.str.contains('Chapter')]
        multi_code_values = multi_code_values.loc[~multi_code_values.Value.str.contains('Block')]

        # if there are any potential ordinal recodes
        if any(ordinal_codes.dataCode == multi_code):
            if multi_code in ordinal_codes.dataCode:
                # if there are recodings
                ordinal_info = ordinal_codes.loc[ordinal_codes.dataCode == multi_code]
                if ordinal_info.reassignments.notnull().values[0]:
                    # first recode any null or edge values
                    multi_code_values = recode_values(df=multi_code_values, reassign=ordinal_info.reassignments)

        # have each level of the cat-multiple as its own row
        multi_out = pd.DataFrame({'fieldid': multi_field, 'array': multi_code_values.Value, 'array_meaning': multi_code_values['Meaning']})
        multi_out = multi_out.merge(ukb_meta, on='fieldid')
        for vis in multi_visits:
            vis_multi_out = copy.deepcopy(multi_out)
            vis_multi_out['visit'] = visit
            list_of_multi_metadata.append(vis_multi_out)

    # merge expanded cat-multi with base metadata df
    multi_metadata_df = pd.concat(list_of_multi_metadata)
    return multi_metadata_df






def combine_phewas_results(config_json, args):
    # combined PHESANT regression results

    # binary logistic
    regr_files = glob.glob(os.path.join(config_json['base_dir'], 'data/ukb/phesant', 'phewas_visit*-results-logistic-binary-*-30.txt'))
    binary_df_list = []
    for stats_file in regr_files:
        binary_df_list.append(read_phewas_stats(stats_file))
    bin_df = pd.concat(binary_df_list, 0)

    # linear
    linear_files = glob.glob(os.path.join(config_json['base_dir'], 'data/ukb/phesant', 'phewas_visit*-results-linear-*-30.txt'))
    linear_df_list = []
    for stats_file in linear_files:
        linear_df_list.append(read_phewas_stats(stats_file))
    linear_df = pd.concat(linear_df_list, 0)

    # ordered logistic
    ordered_logistic_files = glob.glob(os.path.join(config_json['base_dir'], 'data/ukb/phesant', 'phewas_visit*-ordered-logistic-*-30.txt'))
    ordered_df_list = []
    for stats_file in ordered_logistic_files:
        ordered_df_list.append(read_phewas_stats(stats_file))
    ordered_df = pd.concat(ordered_df_list, 0)

    # ordered logistic
    multi_logistic_files = glob.glob(os.path.join(config_json['base_dir'], 'data/ukb/phesant', 'phewas_visit*-multinomial-logistic-*-30.txt'))
    multi_df_list = []
    for stats_file in multi_logistic_files:
        multi_df_list.append(read_phewas_stats(stats_file))
    multi_df = pd.concat(multi_df_list, 0)

    phesant_df = pd.concat([bin_df, linear_df, multi_df, ordered_df])

    # create clean "field" and "array" columns
    phesant_df['tmp_clean'] = phesant_df['id_vis'].replace('x', '')
    phesant_df['fieldid'] = [x.split('_')[0] for x in phesant_df['tmp_clean']]
    array = []
    for id in phesant_df['tmp_clean']:
        id_splits = id.split('_')
        if '#' in id:
            array = array + [id.split('#')[-1]]
        elif len(id_splits) == 2:
            array = array + ['0']
        elif len(id_splits) == 3:
            array = array + [id.split('_')[-1]]
        else:
            array = array + ['0']

    assert phesant_df.shape[0] == len(array)
    phesant_df['array'] = array
    if 'tmp_clean' in phesant_df.columns:
        del phesant_df['tmp_clean']

    return phesant_df


def read_phewas_stats(stats_file):
    visit = stats_file.split('visit')[-1].split('-')[0]
    df = pd.read_csv(stats_file, header=None)
    if df.shape[1] != 1:
        if '-logistic-binary-' in stats_file:
            df.columns = ['id', 'id_vis', 'type', 'case_n', 'beta', 'se_lower', 'se_upper', 'pvalue']
            df.insert(0, 'visit', visit)
            df['odds_ratio'] = np.exp(df['beta'])
            df['n'] = [int(x.split('(')[-1].replace(')', '')) for x in df['case_n']]
            return df

        elif '-linear-' in stats_file or '-ordered-logistic-' in stats_file:
            df.columns = ['id', 'id_vis', 'type', 'n', 'beta', 'se_lower', 'se_upper', 'pvalue']
            df.insert(0, 'visit', visit)
            return df

        elif '-multinomial-logistic-' in stats_file:
            df.columns = ['id', 'id_vis', 'type', 'n', 'beta', 'se_lower', 'se_upper', 'pvalue']
            df.insert(0, 'visit', visit)
            return df


if __name__ == "__main__":
    main()

