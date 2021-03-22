library(tidyverse)
library(rjson)

# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging'
config = fromJSON(file = "input.json")


ref_dir  = paste0(base_dir, '/ref_files/ukb')
out_dir  = paste0(base_dir, '/output/sqlite_files')

# read variable metadata files
showcase_df  = read.csv('/gpfs/milgram/project/holmes/kma52/buckner_aging/ref_files/showcase.csv')
codings_df   = read.csv('/gpfs/milgram/project/holmes/kma52/buckner_aging/ref_files/codings.csv')
#readcodes_df = read.csv(paste0(ref_dir, '/readcodes.csv'))
outcome_df   = read_delim('/gpfs/milgram/project/holmes/kma52/buckner_aging/external/PHESANT/variable-info/outcome-info.tsv', delim='\t')


# pull out relevant fields in each metadata file
df = outcome_df[c('FieldID', 'Category', 'Path', 'Field', 'ValueType', 'DATA_CODING', 'EXCLUDED', 'CAT_MULT_INDICATOR_FIELDS', 'CAT_SINGLE_TO_CAT_MULT', 'DATE_CONVERT')]
showcaseMerge_df = showcase_df[c('FieldID', 'Participants', 'Items', 'Units', 'Strata', 'Instances', 'Array', 'Notes', 'Link')]

# merge them
meta_df = merge(x=df, y=showcaseMerge_df, by='FieldID')

# uniform naming convention
colnames(meta_df) = gsub('EXCLUDED', 'Excluded', colnames(meta_df))
colnames(meta_df) = gsub('DATA_CODING', 'DataCoding', colnames(meta_df))
colnames(meta_df) = gsub('CAT_MULT_INDICATOR_FIELDS', 'CatMultIndicatorFields', colnames(meta_df))
colnames(meta_df) = gsub('DATE_CONVERT', 'DateConvert', colnames(meta_df))
colnames(meta_df) = gsub('CAT_SINGLE_TO_CAT_MULT', 'CatSingleToCatMult', colnames(meta_df))



# Identify the most granular phenotype category
meta_df$phenoCategory = unlist(lapply(meta_df$Path, function(x) strsplit(x, ' > ')[[1]][[ length(strsplit(x, ' > ')[[1]]) ]]))


# Key/Value Pairs
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
ukbVarName = rownames(covar_map_df)
colnames(covar_map_df) = 'ukbField'
covar_map_df[['ukbVarName']] = ukbVarName


# Default Covarate Mapping
covar_dict = list()
covar_dict[['centerSex']] = c(covar_map[['sex']], covar_map[['center']])
covar_dict[['centerSexHeadXYZ']]    = c(covar_dict[['centerSex']], covar_map[['headPosX']], covar_map[['headPosY']], covar_map[['headPosZ']])
covar_dict[['centerSexHeadXYZ_t1NoETIV']] = c(covar_dict[['centerSexHeadXYZ']], covar_map[['t1CNR']], covar_map[['t1SNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV']]       = c(covar_dict[['centerSexHeadXYZ']], covar_map[['eTIV']], covar_map[['t1CNR']], covar_map[['t1SNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV_rfMRI']]    = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['rfmriMotion']], covar_map[['rfmriInvSNR']])
covar_dict[['centerSexHeadXYZ_t1eTIV_dMRI']]     = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['dmriOutliers']])
covar_dict[['centerSexHeadXYZ_t1eTIV_task']]     = c(covar_dict[['centerSexHeadXYZ_t1eTIV']], covar_map[['tfmriMotion']])


# start to match covariate types to imaging fields
meta_df[['covarType']] = 'centerSex'
meta_df[['covarType']][grepl('Diffusion brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV_dMRI'
meta_df[['covarType']][grepl('Task functional brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV_task'
meta_df[['covarType']][grepl('T2-weighted brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV'
meta_df[['covarType']][grepl('Susceptibility weighted brain MRI', meta_df$Path)] = 'centerSexHeadXYZ_t1eTIV'

t1_matches = grep('T1 structural brain MRI', meta_df$Path)
vol_matches = grep('Volume of ', meta_df$Field)
t1_vol_matches = intersect(t1_matches, vol_matches)
meta_df[['covarType']][t1_vol_matches] = 'centerSexHeadXYZ_t1eTIV'

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


# write output
write.csv(covariate_df, paste0(out_dir, '/covariateTable.csv'), row.names=F, quote=T)
write.csv(meta_df, paste0(out_dir, '/ukbMetaData.csv'), row.names=F, quote=T)













