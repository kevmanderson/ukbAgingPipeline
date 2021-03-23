library(RSQLite)
library(reshape2)
library(tidyverse)
library(rjson)
library(optparse)
options(stringsAsFactors=F)

option_list = list(
  make_option(c("-c", "--config"), type="character", default=NULL, dest='config',
              help="full path to config file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# process config file
config   = fromJSON(file = opt$config)
base_dir = config[[1]]$base_dir
showcase_file = config[[1]]$showcase_file
outcome_file  = config[[1]]$outcome_file
codings_file  = config[[1]]$codings_file
ordinal_file  = config[[1]]$ordinal_file


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

# sql file to create
sqlite3_file = paste0(out_dir, '/ukb_sqlite3.db')
sqlite3_file = paste0(out_dir, '/ukb_sqlite3_100k.db')

# files created by ukb_conv
script_path = gsub('.enc_ukb', '.r', config[[1]]$ukb_enc)
tab_file    = gsub('.enc_ukb', '.tab', config[[1]]$ukb_enc)
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
con    = file(tab_file, "r")
header = as.character(read.table(con, header=F, nrows=1, sep="\t"))

# create metadata df that matches columns of ukb df
ukbFields  = unlist(lapply(header, function(x) strsplit(x,'[.]')[[1]][[2]]))
ukbFieldDF = data.frame(FieldID=ukbFields, ColName=header)
ukbFieldDF = merge(x=ukbFieldDF, y=ukbMetaData, by='FieldID')

# fields that don't require recoding
nonMultFields = which(is.na(ukbFieldDF$CatMultIndicatorFields))

if (file.exists(sqlite3_file)){
  file.remove(sqlite3_file)
}
conn = dbConnect(RSQLite::SQLite(), sqlite3_file)

# chunk through and reformat data 10k at a time
chunk = 1
chunk_size = 10000
while (con){
  # while (chunk <= 10){
  write(chunk, '')
  # read UKB data chunk
  bd = read.table(con, header=F, nrows=chunk_size, sep="\t")
  colnames(bd) = header

  # call the script that was created above. Convert numeric to strings
  source(script_outPath)
  chunk = chunk + 1

  # write data to different SQL tables, split by data type
  dtypePairs = list(c('int','Integer'), c('real', 'Continuous'), c('str','Categorical single'), c('str','Categorical multiple'), c('date', 'Date'))
  pair = dtypePairs[[1]]
  for (pair in dtypePairs){
    write(pair,'')
    dtype   = pair[[2]]
    sqlType = pair[[1]]

    # find columns for this data type (make sure they dont have multi-array codings)
    dtype_idxs    = intersect(which(ukbFieldDF$ValueType == dtype), nonMultFields)
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

dbDisconnect(conn)

close(conn)
close(con)








# write data to different SQL tables
dtypePairs = list(c('int','Integer'), c('real', 'Continuous'), c('str','Categorical single'), c('str','Categorical multiple'), c('Date', 'date'))
pair = dtypePairs[[1]]
for (pair in dtypePairs){
  dtype   = pair[[2]]
  sqlType = pair[[1]]

  # find columns for this data type
  dtype_idxs    = intersect(which(ukbFieldDF$ValueType == dtype), nonMultFields)
  dtypeFieldIds = c('f.eid', ukbFieldDF$ColName[dtype_idxs])
  dtype_df      = ukb_df[as.character(dtypeFieldIds)]

  # melt data to long before writing to SQL
  dtype_melted = reshape2::melt(dtype_df, id.vars=c('f.eid'), na.rm=T)
  colnames(dtype_melted) = c('eid', 'variable', 'value')

  # split column name into
  dtype_melted$variable = as.character(dtype_melted$variable)
  var_splits   = lapply(dtype_melted$variable, function(x) strsplit(x,'[.]')[[1]])
  dtype_melted$field = unlist(lapply(var_splits, function(x) x[[2]]))
  dtype_melted$time  = unlist(lapply(var_splits, function(x) x[[3]]))
  dtype_melted$array = unlist(lapply(var_splits, function(x) x[[4]]))
  dtype_melted$variable = gsub('[.]','_', gsub('f.', '', dtype_melted$variable))
  #dtype_melted$variable = NULL

  dbWriteTable(conn, sqlType, dtype_melted, append = TRUE)
}


dtype_melted = reshape2::melt(dtype_df, id.vars=c('f.eid'), na.rm=T)
dtype_melted$variable = as.character(dtype_melted$variable)
var_splits   = lapply(dtype_melted$variable, function(x) strsplit(x,'[.]')[[1]])
dtype_melted$field = unlist(lapply(var_splits, function(x) x[[2]]))
dtype_melted$time  = unlist(lapply(var_splits, function(x) x[[3]]))
dtype_melted$array = unlist(lapply(var_splits, function(x) x[[4]]))
dtype_melted$variable = NULL


colnames(dtype_melted) = c('eid','variable','value')
dbWriteTable(conn, 'int', dtype_melted, append = TRUE)

# write metadata to sql
conn = dbConnect(RSQLite::SQLite(), sqlite3_file)
dbWriteTable(conn, 'metadata', ukbFieldDF, replace = FALSE)


dbDisconnect(conn)



library(data.table)
long <- melt(setDT(wide), id.vars = c("Code","Country"), variable.name = "year")




