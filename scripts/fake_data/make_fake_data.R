
library(RSQLite)
library(tidyverse)
library(rjson)
library(optparse)
library(dplyr)
options('stringsAsFactors'=F)

option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, dest='config',
                help="full path to config file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# repo_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline'
# base_dir = "/gpfs/milgram/project/holmes/kma52/buckner_aging"
# process config file
config   = fromJSON(file = opt$config)
base_dir = config[[1]]$base_dir
repo_dir = config[[1]]$repo_dir
showcase_file = paste0(repo_dir, '/ref_files/showcase.csv')
outcome_file  = paste0(repo_dir, '/ref_files/outcome-info.tsv')
codings_file  = paste0(repo_dir, '/ref_files/data-coding-ordinal-info.txt')

# output directory for synthetic data
# --------
fake_ref_dir = paste0(base_dir, '/data/fake_data')


# read variable metadata files
showcase_df  = read.csv(showcase_file)
codings_df   = read.csv(codings_file)
outcome_df   = read_delim(outcome_file, delim='\t')

# pull out relevant fields in each metadata file
df = outcome_df[c('FieldID', 'Category', 'Path', 'Field', 'ValueType', 'DATA_CODING', 'EXCLUDED', 'CAT_MULT_INDICATOR_FIELDS', 'CAT_SINGLE_TO_CAT_MULT', 'DATE_CONVERT')]
showcaseMerge_df = showcase_df[c('FieldID', 'Participants', 'Items', 'Units', 'Strata', 'Instances', 'Array', 'Notes', 'Link')]

# merge them
meta_df = merge(x=df, y=showcaseMerge_df, by='FieldID')

# Change the name formatting for some PHESANT metadata fields
colnames(meta_df) = gsub('EXCLUDED', 'Excluded', colnames(meta_df))
colnames(meta_df) = gsub('DATA_CODING', 'DataCoding', colnames(meta_df))
colnames(meta_df) = gsub('CAT_MULT_INDICATOR_FIELDS', 'CatMultIndicatorFields', colnames(meta_df))
colnames(meta_df) = gsub('DATE_CONVERT', 'DateConvert', colnames(meta_df))
colnames(meta_df) = gsub('CAT_SINGLE_TO_CAT_MULT', 'CatSingleToCatMult', colnames(meta_df))

# Identify the most high-level phenotype category
meta_df$phenoCategory = unlist(lapply(meta_df$Path, function(x) strsplit(x, ' > ')[[1]][[ length(strsplit(x, ' > ')[[1]]) ]]))

# Key/Value Pairs. Link field descriptors to their UKB showcase ids
covar_map = list()
covar_map[['sex']] = '22001'
covar_map[['center']] = '54'
covar_map[['t1CNR']] = '25735'
covar_map[['t1SNR']] = '25734'
covar_map[['headPosX']] = '25756'
covar_map[['headPosY']] = '25757'
covar_map[['headPosZ']] = '25758'
covar_map[['eTIV']] = '26521'
covar_map[['rfmriMotion']] = '25741'
covar_map[['rfmriInvSNR']] = '25744'
covar_map[['tfmriMotion']] = '25742'
covar_map[['dmriOutliers']] = '25746'
covar_map_df = as.data.frame(t(as.data.frame(covar_map)))

# format the resulting dataframe (ukbField, ukbVarName are the two columns)
ukbVarName = rownames(covar_map_df)
colnames(covar_map_df) = 'ukbField'
covar_map_df[['ukbVarName']] = ukbVarName
# head(covar_map_df)

# Default Covariate Mapping
# different imaging categories will have different initial covariate sets
covar_dict = list()
covar_dict[['centerSex']] = c(covar_map[['sex']], covar_map[['center']])
covar_dict[['centerSexHeadXYZ']]    = c(covar_dict[['centerSex']], covar_map[['headPosX']], covar_map[['headPosY']], covar_map[['headPosZ']])
covar_dict[['centerSexHeadXYZ_t1NoETIV']] = c(covar_dict[['centerSexHeadXYZ']], covar_map[['t1CNR']], covar_map[['t1SNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV']]       = c(covar_dict[['centerSexHeadXYZ']], covar_map[['eTIV']], covar_map[['t1CNR']], covar_map[['t1SNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV_rfMRI']] = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['rfmriMotion']], covar_map[['rfmriInvSNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV_dMRI']]  = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['dmriOutliers']])
covar_dict[['centerSexHeadXYZ_t1eTIV_task']]  = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['tfmriMotion']])


# start to match covariate types to specific imaging fields

# everything starts as center/sex
meta_df[['covarType']] = 'centerSex'

# diffusions/task/t2/swi fields all control for t1-eTIV
meta_df[['covarType']][grepl('Diffusion brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV_dMRI'
meta_df[['covarType']][grepl('Task functional brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV_task'
meta_df[['covarType']][grepl('T2-weighted brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV'
meta_df[['covarType']][grepl('Susceptibility weighted brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV'

# t1/freesurfer structurals require a bit more fine-grained matching to covariates
t1_matches = grep('T1 structural brain MRI', meta_df$Path)

# volumetric measures control for t1/eTIV
vol_matches = grep('Volume of ', meta_df$Field)
t1_vol_matches = intersect(t1_matches, vol_matches)
meta_df[['covarType']][t1_vol_matches] = 'centerSexHeadXYZ_t1eTIV'

# thickness/area/grey-white all don't control for eTIV
nonVol_matches = grep('Mean thickness of |Area of |Grey-white contrast |Mean intensity of ', meta_df$Field)
t1_nonVol_matches = intersect(t1_matches, nonVol_matches)
meta_df[['covarType']][t1_nonVol_matches] = 'centerSexHeadXYZ_t1NoETIV'


# create a table with covariate information
covariate_df = NULL
for (name in names(covar_dict)){
  write(name,'')
  addDat = covar_map_df[covar_map_df$ukbField %in% covar_dict[[name]],]
  addDat$covarType = name
  covariate_df = rbind(covariate_df, addDat)
}
rownames(covariate_df) = NULL


# start to make the fake dataframe
nsubs = 10000
synethic_df = data.frame(f.eid=1:nsubs)

# make fake age data
synethic_df['f.21003.0.0'] = rnorm(n=nsubs, mean=65, sd=5)
synethic_df['f.21003.1.0'] = rnorm(n=nsubs, mean=65, sd=5)
synethic_df['f.21003.2.0'] = rnorm(n=nsubs, mean=65, sd=5)
synethic_df['f.21003.3.0'] = rnorm(n=nsubs, mean=65, sd=5)

# genetic sex
synethic_df['f.22001.0.0'] = ifelse(sample(c(0,1), replace=TRUE, size=nsubs) == 1, 'Male', 'Female')

# UKB assessment center
center_arr = sample(c(0,1,2), replace=TRUE, size=nsubs)
synethic_df['f.54.2.0'] = recode(center_arr, `0`='Cheadle (imaging)', `1`='Newcastle (imaging)', `2`='Reading (imaging)')
center_arr = sample(c(0,1,2), replace=TRUE, size=nsubs)
synethic_df['f.54.0.0'] = recode(center_arr, `0`='Cheadle', `1`='Newcastle', `2`='Reading')

# structural MRI covariates
synethic_df['f.25735.2.0'] = rnorm(n=nsubs, mean=80, sd=5)
synethic_df['f.25734.2.0'] = rnorm(n=nsubs, mean=30, sd=2)
synethic_df['f.25756.2.0'] = rnorm(n=nsubs, mean=10, sd=12)
synethic_df['f.25757.2.0'] = rnorm(n=nsubs, mean=100, sd=15)
synethic_df['f.25758.2.0'] = rnorm(n=nsubs, mean=40, sd=4)
synethic_df['f.26521.2.0'] = rnorm(n=nsubs, mean=1000, sd=4)

# functional MRI covariates
synethic_df['f.25741.2.0'] = rnorm(n=nsubs, mean=.3, sd=.01)
synethic_df['f.25742.2.0'] = rnorm(n=nsubs, mean=.28, sd=.02)
synethic_df['f.25744.2.0'] = rnorm(n=nsubs, mean=1, sd=.3)

# genetic principle components
for (pc in 0:39){
  synethic_df[paste0('f.22009.0.', pc)] = rnorm(n=nsubs, mean=0, sd=.5)
}

# make fake freesurfer data
freesurfer_cols = meta_df[meta_df$phenoCategory == 'Freesurfer ASEG',]$FieldID
for (field_id in freesurfer_cols){
  synethic_df[paste0('f.', field_id, '.0.0')] = rnorm(n=nsubs, mean=sample(2000:4000, 1), sd=sample(20:40, 1))
}

header = colnames(synethic_df)
ukbMetaData = meta_df

# create metadata df that matches columns of ukb df
ukbFields  = unlist(lapply(header, function(x) strsplit(x,'[.]')[[1]][[2]]))
ukbFieldDF = data.frame(FieldID=ukbFields, ColName=header)
ukbFieldDF = merge(x=ukbFieldDF, y=ukbMetaData, by='FieldID')

ukbFieldDF$visit = unlist(lapply(as.character(ukbFieldDF$ColName), function(x) strsplit(x,'[.]')[[1]][[3]]))
ukbFieldDF$visit = as.numeric(ukbFieldDF$visit)


# create sqlite file
fake_sql_file = paste0(fake_data_dir, 'ukb_synthetic_sqlite3.db')
if (file.exists(fake_sql_file)){
  file.remove(fake_sql_file)
}
conn = dbConnect(RSQLite::SQLite(), fake_sql_file)

# write data to different SQL tables, split by data type
dtypePairs = list(c('int','Integer'), c('real', 'Continuous'), c('str','Categorical single'), c('str','Categorical multiple'), c('date', 'Date'))
pair = dtypePairs[[2]]
for (pair in dtypePairs){
  write(pair,'')
  dtype   = pair[[2]]
  sqlType = pair[[1]]

  # find columns for this data type (make sure they dont have multi-array codings)
  dtype_idxs    = which(ukbFieldDF$ValueType == dtype)
  dtypeFieldIds = c('f.eid', ukbFieldDF$ColName[dtype_idxs])
  if (length(dtypeFieldIds) == 1){
    next
  }
  dtype_df      = synethic_df[as.character(dtypeFieldIds)]

  # melt data to long before writing to SQL
  dtype_melted = reshape2::melt(dtype_df, id.vars=c('f.eid'), na.rm=T)
  colnames(dtype_melted) = c('eid', 'variable', 'value')

  # split column name into
  dtype_melted$variable = as.character(dtype_melted$variable)
  var_splits   = lapply(dtype_melted$variable, function(x) strsplit(x,'[.]')[[1]])
  dtype_melted$field = as.numeric(unlist(lapply(var_splits, function(x) x[[2]])))
  dtype_melted$time  = as.numeric(unlist(lapply(var_splits, function(x) x[[3]])))
  dtype_melted$array = as.numeric(unlist(lapply(var_splits, function(x) x[[4]])))
  dtype_melted$variable = gsub('[.]','_', gsub('f.', '', dtype_melted$variable))
  #dtype_melted$variable = NULL
  #fieldTypes = c('eid'='INTEGER', 'variable'='TEXT', 'value'='DATE', 'field'='INTEGER', 'time'='INTEGER', 'array'='INTEGER')
  #CREATE TABLE `date` ( `eid` INTEGER, `variable` TEXT, `value` REAL, `field` TEXT, `time` TEXT, `array` TEXT )
  dbWriteTable(conn, sqlType, dtype_melted, append=TRUE)
}
# index each of the created sql tables
dbExecute(conn, 'CREATE INDEX str_index ON str (variable, value, eid)')
dbExecute(conn, 'CREATE INDEX int_index ON int (variable, value, eid)')
dbExecute(conn, 'CREATE INDEX real_index ON real (variable, value, eid)')
#dbExecute(conn, 'CREATE INDEX dt_index ON date (variable, value, eid)')


# format column
ukbFieldDF$ColName = gsub('f[.]','',ukbFieldDF$ColName)
ukbFieldDF$ColName = gsub('[.]','_',ukbFieldDF$ColName)

# add study visit
ukbFieldDF$Visit = as.numeric(unlist(lapply(ukbFieldDF$ColName, function(x) strsplit(x,'_')[[1]][[2]])))

# add some fields to metadata
ukbFieldDF$SqlType = NA
ukbFieldDF$tab = NA
dtypePairs = list(c('int','Integer','INTEGER'), c('real', 'Continuous','REAL'), c('str','Categorical single','VARCHAR'), c('str','Categorical multiple','VARCHAR'), c('date', 'Date','NUMERIC'))
for (pair in dtypePairs){
  dtype   = pair[[2]]
  sqlType = pair[[3]]
  tab     = pair[[1]]
  # find columns for this data type
  dtypeColNames = ukbFieldDF$ColName[which(ukbFieldDF$ValueType == dtype)]
  if (length(dtypeColNames) == 0){
    next
  }
  ukbFieldDF$SqlType[ukbFieldDF$ColName %in% dtypeColNames] = sqlType
  ukbFieldDF$tab[ukbFieldDF$ColName %in% dtypeColNames] = tab
}

# commit the other data tables
dbWriteTable(conn, 'metadata', ukbFieldDF, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'codingsDF', outcome_df, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'covariate_table', covariate_df, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'code_ordinal', codings_df, overwrite = TRUE, append=FALSE)

dbDisconnect(conn)



