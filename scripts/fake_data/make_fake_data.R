
library(tidyversese)

repo_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline'
fake_ref_dir = paste0(repo_dir, '/scripts/fake_data')

# read ukb metadata
meta_path = paste0(fake_ref_dir, '/ref/ukbMetaData.csv')
meta_df = read.csv(meta_path)

codings_df = read.csv(paste0(repo_dir, '/ref_files/codings.csv'))

nsubs = 1000

# don't bother with dates
meta_df = meta_df[meta_df$ValueType != 'Date',]

# or complicated fields
meta_df = meta_df[is.na(meta_df$CatMultIndicatorFields),]
meta_df = meta_df[is.na(meta_df$CatSingleToCatMult),]
meta_df = meta_df[!meta_df$Excluded %in% c('YES-PROCESSING', 'YES-POLYMORPHIC', 'YES-CAT-SIN-MUL-VAL', 'YES-SENSITIVE', 'YES-ACE'),]

fake_matrix = matrix(ncol=nrow(meta_df), nrow=nsubs)
for (row in 1:nrow(meta_df)){
  meta_row = meta_df[row,]
  if (meta_row$ValueType == 'Integer'){
    fake_dat = sample(0:100, nsubs, replace=T)
  } else if (meta_row$ValueType == 'Continuous'){
    fake_dat = sample(0:10000, nsubs, replace=T)/100
  } else if (meta_row$ValueType == 'Categorical single' || meta_row$ValueType == 'Categorical Multiple'){
    write(row,'')
    meanings = as.character(codings_df$Meaning[codings_df$Coding == meta_row$DataCoding])
    meanings = as.factor(meanings)
    fake_levels = sample(1:length(unique(meanings)), nsubs, replace=T)
    fake_dat = meanings[fake_levels]
  }
  fake_matrix[,row] = fake_dat
}
fake_df = as.data.frame(fake_matrix)






