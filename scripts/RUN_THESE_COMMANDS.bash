#!/bin/bash


# set up paths
# -------
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
repo_dir=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
data_dir=/gpfs/milgram/project/holmes/kma52/buckner_aging


# 1: Create Directories
# --------
python3 ./scripts/ukb_create_dirs.py \
  --root_dir=${data_dir}


# 2: Decrypt UKB data
# --------
python3 ./scripts/ukb_decrypt.py \
  --enc_key_pair=${data_dir}/data/ukb/raw/ukb40501.enc:${data_dir}/data/ukb/raw/ukb40501.key \
  --out_dir=${data_dir}


# 2: Convert UKB data
# --------
python3 ./scripts/ukb_convert.py \
  --enc_ukb=${data_dir}/data/ukb/raw/ukb40501.enc_ukb \
  --ukb_conv=${data_dir}/data/ukb/external/ukbconv \
  --fields_to_convert=${repo_dir}/ref_files/fields_to_convert.txt \
  --out_dir=${data_dir}/data/ukb/raw


# 3: Make bulk lists
# --------
python3 ./scripts/ukb_make_bulk_lists.py \
  --ukb_conv=${repo_dir}/external/ukb_tools/ukbconv \
  --enc_ukb=${data_dir}/data/ukb/raw/ukb40501.enc_ukb \
  --enc_key=${data_dir}/data/ukb/raw/ukb40501.key \
  --bulk_field='rfmri_full_25:25750'  \
  --bulk_field='rfmri_full_100:25751'  \
  --bulk_field='rfmri_part_25:25752'  \
  --bulk_field='rfmri_part_100:25753'  \
  --bulk_field='rfmri_rsfa_25:25754'  \
  --bulk_field='rfmri_rsfa_100:25755' \
  --out_dir=${data_dir}/data/ukb/bulk \
  --slurm_dir=${data_dir}/slurm


# 4: Download bulk lists
# --------
python ./scripts/ukb_download_bulk_data.py \
    --ukb_conv=${repo_dir}/external/ukb_tools/ukbfetch \
    --enc_key=${data_dir}/data/ukb/raw/ukb40501.key \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755' \
    --out_dir=${data_dir}/data/ukb/bulk \
    --slurm_dir=${data_dir}/slurm


# 5: Bulk download to csv
# --------
python ./scripts/ukb_download_bulk_data.py \
    --ukb_conv=${repo_dir}/external/ukb_tools/ukbfetch \
    --enc_key=${data_dir}/data/ukb/raw/ukb40501.key \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755' \
    --out_dir=${data_dir}/data/ukb/bulk \
    --slurm_dir=${data_dir}/slurm



# 6: PREPARE for phesant
# --------
python3 ./scripts/ukb_phesant_prepare.py \
    --out_dir=${data_dir}/data/ukb/raw \
    --outcome_info=${repo_dir}/ref_files/outcome-info.tsv \
    --field_file=${repo_dir}/ref_files/field.txt \
    --showcase=${repo_dir}/ref_files/Data_Dictionary_Showcase.tsv \
    --phesant_csv=${data_dir}/data/ukb/raw/ukb40501.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25750_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25750_3.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25751_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25751_3.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25752_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25752_3.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25753_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25753_3.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25754_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25754_3.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25755_2.csv \
    --phesant_csv=${data_dir}/data/ukb/raw/bulk_25755_3.csv \
    --phesant_visits='0;1;2;3' \
    --write_metadata=${repo_dir}/ref_files/aging-outcome-info.tsv


# 7: RUN phesant
# --------
# example for visit-0
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline

age_var=x21003_0_0
python3 ./scripts/ukb_phesant_run.py \
  --pheno_file=${data_dir}/data/ukb/raw/phesant_visit0.csv \
  --confounder_file=${data_dir}/data/ukb/raw/phesant_covar_vis0.csv \
  --metadata_file=${repo_dir}/ref_files/aging-outcome-info.tsv \
  --datacoding_file=${repo_dir}/ref_files/data-coding-ordinal-info.txt \
  --phesant_dir=${repo_dir}/external/PHESANT/WAS \
  --age_var=${age_var} \
  --visit=0 \
  --nparts=30 \
  --out_dir=${data_dir}/data/ukb/phesant \
  --write_scripts_dir=${data_dir}/slurm

  

# 7: Preprocess imputed SNP data
# --------
snp_dir=/gpfs/milgram/data/UKB/ukb_snp
plink2=${data_dir}/external/plink2
python ./scripts/ukb_snp_preprocess.py \
    --ukb_csv=${data_dir}/data/ukb/raw/phesant_visit0.csv \
    --bgens=${snp_dir}/ukb_imp_chr*_v3.bgen \
    --samples=${snp_dir}/ukb25163_imp_chr*_v3_s487324.sample \
    --info_files=${snp_dir}/ukb_mfi_chr*_v3.txt \
    --plink2={plink2} \
    --maf=0.01 \
    --hwe=1e-6 \
    --mind=0.1 \
    --geno=0.1 \
    --info=0.8 \
    --out_dir=${data_dir}/data/ukb/genetic/imputed_processed \
    --slurm_dir=${data_dir}/slurm \
    --slurm


# 8: Polygenic Score Calculation
# --------
plink2=/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2
python ./scripts/ukb_polygenic \

# 8: Polygenic Score Calculation
# --------
plink2=/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2
python ./scripts/ukb_polygenic \








