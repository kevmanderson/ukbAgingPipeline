library(RSQLite)
library(reshape2)
library(tidyverse)
library(rjson)
library(optparse)
options(stringsAsFactors=F)


recode_mult_fields = function(bd, convert_field){
  grep_me = paste0('f.',convert_field, '.')
  convert_data = bd[grep(grep_me, colnames(bd))]
  uniq_visits = unique(unlist(lapply(colnames(convert_data), function(x) strsplit(x, '[.]')[[1]][[3]])))
  uniq_vis = uniq_visits[1]

  new_df = NULL
  for (uniq_vis in uniq_visits){
    vis_grep_me = paste0('f.',convert_field, '.', uniq_vis)
    vis_convert_data = convert_data[grep(grep_me, colnames(convert_data))]

    # keep NA values as NA
    all_na_arr = apply(is.na(vis_convert_data), 1, all)

    uniq_values = unique(unlist(apply(vis_convert_data, 2, unique)))
    uniq_values = uniq_values[!is.na(uniq_values)]

    if (length(uniq_values) == 0){
      next
    } else {
      vis_new_df = NULL
      for (uniq_value in uniq_values){
        new_col_name = paste0(vis_grep_me, '#', uniq_value)
        binary_indicator = rowSums(vis_convert_data == uniq_value, na.rm=T)
        # set NAs
        binary_indicator[all_na_arr] = NA
        vis_new_df[[new_col_name]] = binary_indicator
      }
      vis_new_df = data.frame(vis_new_df)
      new_df = c(new_df, vis_new_df)
    }
  }
  if (length(new_df) == 0){
    filter_bd = bd[,!grepl(grep_me, colnames(bd))]
    return (filter_bd)
  } else {
    new_df = do.call('cbind', new_df)
    new_df = as.data.frame(new_df)

    filter_bd = bd[,!grepl(grep_me, colnames(bd))]
    filter_bd = cbind(filter_bd, new_df)

    return(filter_bd)
  }
}



option_list = list(
  make_option(c("-c", "--config"), type="character", default=NULL, dest='config',
              help="full path to config file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt = NULL
opt$config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'

# process config file
config   = fromJSON(file = opt$config)
base_dir = config[[1]]$base_dir
showcase_file = '/ref_files/showcase.csv' #config[[1]]$showcase_file
outcome_file  = '/ref_files/outcome-info.tsv' #config[[1]]$outcome_file
codings_file  = '/ref_files/codings.csv' #config[[1]]$codings_file
ordinal_file  = '/ref_files/data-coding-ordinal-info.txt' #config[[1]]$ordinal_file


# set up directories
#base_dir = '/ncf/sba01/brain_aging'
#ref_dir  = paste0(base_dir, '/ref_files/ukb')
out_dir  = paste0(base_dir, '/output/sqlite_files')

# Read MetaData info
codingsDF   = read.csv(codings_file)
ukbMetaData = read.csv(paste0(out_dir, '/ukbMetaData.csv'))
ukbMetaData$FieldID = as.character(ukbMetaData$FieldID)

# previosuly created covariate table
covariateTable = read.csv(paste0(out_dir, '/covariateTable.csv'))

# ordinal data field information
ordinalDF = read.csv(ordinal_file)

# sql file to create
sqlite3_file = paste0(out_dir, '/ukb_sqlite3.db')
sqlite3_file = paste0(out_dir, '/ukb_sqlite3_100k.db')

# files created by ukb_conv
script_path = gsub('.enc', '.r', config[[1]]$ukb_enc)
tab_file    = gsub('.enc', '.tab', config[[1]]$ukb_enc)
#script_path = '/ncf/sba01/brain_aging/data/ukb/raw/ukb45156.r'
#tab_file = '/ncf/sba01/brain_aging/data/ukb/raw/ukb45156.tab'

# re-write r script to keep the label recodings only.
# removes lines that actually reads the data, we instead chunk through the raw data below.
script_outPath = gsub('[.]r', '_recode.r', script_path)
res = readLines(script_path)
writeThese = NULL
for (line in res){
  # don't read the one line that tries to load all data at once
  if ( !grepl('bd <- ', line) ) {
    writeThese = c(writeThese, line)
  }
}
# write new script
fileConn = file(script_outPath)
writeLines(writeThese, fileConn)
close(fileConn)
# --------

# read header info
#con    = file(tab_file, "r")
#header = as.character(read.table(tab_file, header=F, nrows=1, sep="\t"))

h = read.table(tab_file, header=F, nrows=1, sep="\t")
header = as.character(t(h))

# create metadata df that matches columns of ukb df
ukbFields  = unlist(lapply(header, function(x) strsplit(x,'[.]')[[1]][[2]]))
ukbFieldDF = data.frame(FieldID=ukbFields, ColName=header)
ukbFieldDF = merge(x=ukbFieldDF, y=ukbMetaData, by='FieldID')

# clear out some fields we don't want to examine
ukbFieldDF = ukbFieldDF[!grepl('YES-SUPERSEDED', ukbFieldDF$Excluded),]
ukbFieldDF = ukbFieldDF[!grepl('YES-DATE', ukbFieldDF$Excluded),]
ukbFieldDF = ukbFieldDF[!grepl('YES-POLYMORPHIC', ukbFieldDF$Excluded),]
ukbFieldDF = ukbFieldDF[!grepl('YES-PROCESSING', ukbFieldDF$Excluded),]

ukbFieldDF$visit = unlist(lapply(as.character(ukbFieldDF$ColName), function(x) strsplit(x,'[.]')[[1]][[3]]))
ukbFieldDF$visit = as.numeric(ukbFieldDF$visit)

# medication codes are super messy, manually exclude for now
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '20003'),]
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '20199'),]

# Edge case matrix completion variable
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '6332'),]

# cancer histology
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '40011'),]

# sensitive and not useful
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '41216'),]
ukbFieldDF = ukbFieldDF[which(ukbFieldDF$FieldID != '41215'),]


#ukbFieldDF$ColName = lapply(ukbFieldDF$ColName, function(x) str_replace_all(x, "[[:punct:]]| ", '.'))

x = data.frame(t(matrix(NA, length(ukbFieldDF$ColName))))
colnames(x) = ukbFieldDF$ColName
# format multiple-categorical outcomes
cat_mult_df = ukbFieldDF[which(ukbFieldDF$ValueType == 'Categorical multiple'),]
mult_field_ids = as.character(unique(cat_mult_df$FieldID))
for (field_id in mult_field_ids){
  # unformatted ukb metadata
  mult_fields = ukbFieldDF[which(ukbFieldDF$FieldID == field_id),]

  # data coding info required for creatino of new ukb metadata entries
  cur_codes = codingsDF[codingsDF$Coding %in% mult_fields$DataCoding[1],]

  # get rid of negative non answer
  cur_codes$Value = as.numeric(as.character(cur_codes$Value))
  cur_codes_nonan = cur_codes[cur_codes$Value > 0,]

  # replace the "colname"" separately for each visit, since there are some visit-specific fields taht we want to retain.
  uniq_visits = unique(mult_fields$visit)
  for (vis in uniq_visits){
    write(vis,'')
    template_row  = mult_fields[which(mult_fields$visit == vis)[1],]
    new_name_base = paste0('f.', field_id, '.', as.character(vis))
    new_name_arr  = paste0(new_name_base, '#', as.character(cur_codes_nonan$Meaning))
    new_field_df = NULL
    for (n in new_name_arr){
      print(n)
      template_row$ColName = n
      #new_field_list = c(new_field_list, template_row)
      new_field_df = rbind(new_field_df, template_row)
    }
    # get rid of old rows
    ukbFieldDF = ukbFieldDF[!grepl(field_id, ukbFieldDF$FieldID),]
    # add new rows
    ukbFieldDF = rbind(ukbFieldDF, new_field_df)
  }
}

x = NULL
for (n in ukbFieldDF$ColName){
  x[[n]] = c(0,1,0)
}
{}

# 22599
# 	22617
# employ
#df = ukbFieldDF[ukbFieldDF$FieldID == '22617',]


# fields that don't require recoding
nonMultFields = which(is.na(ukbFieldDF$CatMultIndicatorFields))
convertMultFields = unique(ukbFieldDF$FieldID[!is.na(ukbFieldDF$CatMultIndicatorFields)])

if (file.exists(sqlite3_file)){
  file.remove(sqlite3_file)
}
conn = dbConnect(RSQLite::SQLite(), sqlite3_file)

# read header info
con    = file(tab_file, "r")
header = as.character(read.table(tab_file, header=F, nrows=1, sep="\t"))

# chunk through and reformat data 10k at a time
chunk = 1
chunk_size = 10000
while (con){

  # while (chunk <= 10){
  write(chunk, '')
  # read UKB data chunk
  if (chunk == 1){
    bd = read.table(con, header=T, nrows=chunk_size, sep="\t")
  } else {
    bd = read.table(con, header=F, nrows=chunk_size, sep="\t")
  }
  colnames(bd) = header
  # call the script that was created above. Convert numeric to strings
  source(script_outPath)
  chunk = chunk + 1

  # convert
  convert_field = convertMultFields[1]
  bd_save = bd
  for (convert_field in convertMultFields){
    write(convert_field,'')
    bd = recode_mult_fields(bd=bd, convert_field=convert_field)
  }

  # write data to different SQL tables, split by data type
  dtypePairs = list(c('int','Integer'), c('real', 'Continuous'), c('str','Categorical single'), c('str','Categorical multiple'), c('date', 'Date'))
  pair = dtypePairs[[3]]
  for (pair in dtypePairs){
    write(pair,'')
    dtype   = pair[[2]]
    sqlType = pair[[1]]

    # find columns for this data type (make sure they dont have multi-array codings)
    #dtype_idxs    = intersect(which(ukbFieldDF$ValueType == dtype), nonMultFields)
    dtype_idxs    = which(ukbFieldDF$ValueType == dtype)
    dtypeFieldIds = c('f.eid', ukbFieldDF$ColName[dtype_idxs])
    dtype_df      = bd[as.character(dtypeFieldIds)]

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
}
# index each of the created sql tables
dbExecute(conn, 'CREATE INDEX str_index ON str (variable, value, eid)')
dbExecute(conn, 'CREATE INDEX int_index ON int (variable, value, eid)')
dbExecute(conn, 'CREATE INDEX real_index ON real (variable, value, eid)')
dbExecute(conn, 'CREATE INDEX dt_index ON date (variable, value, eid)')


# TODO: ADD SQL DATATYPE

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
  dtypeColNames = ukbFieldDF$ColName[intersect(which(ukbFieldDF$ValueType == dtype), nonMultFields)]
  ukbFieldDF$SqlType[ukbFieldDF$ColName %in% dtypeColNames] = sqlType
  ukbFieldDF$tab[ukbFieldDF$ColName %in% dtypeColNames] = tab
}

# commit the other data tables
dbWriteTable(conn, 'metadata', ukbFieldDF, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'codingsDF', ukbFieldDF, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'covariate_table', covariateTable, overwrite = TRUE, append=FALSE)
dbWriteTable(conn, 'code_ordinal', ordinalDF, overwrite = TRUE, append=FALSE)

dbDisconnect(conn)
close(conn)
close(con)





