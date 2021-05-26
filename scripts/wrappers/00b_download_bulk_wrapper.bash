#!/bin/bash


cd /ncf/sba01/ukbAgingPipeline

python ./scripts/main.py \
  --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
  --stage='decrypt'

python ./scripts/main.py \
  --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
  --stage='convert'

python ./scripts/main.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
    --stage='make_bulk_list' \
    --bulk-field='rfmri_full_25:25750'  \
    --use-enc-ukb='ukb45156.enc_ukb'

python ./scripts/main.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
    --stage='make_bulk_list' \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755' \
    --use-enc-ukb='ukb45156.enc_ukb'

python ./scripts/main.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
    --stage='make_bulk_list' \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755' \
    --use-enc-ukb='ukb45156.enc_ukb'

cd /ncf/sba01/ukbAgingPipeline
python ./scripts/main.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
    --stage='download_bulk' \
    --bulk-field='mri_t2_nii:20253'

cd /ncf/sba01/ukbAgingPipeline
python ./scripts/main.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json \
    --stage='download_bulk' \
    --bulk-field='mri_rest_nii:20227' \
     --bulk-field='mri_rest_dicom:20225'  \
     --bulk-field='mri_rest_nii:20227'  \
     --bulk-field='mri_t1_nii:20252'  \
     --bulk-field='mri_t2_nii:20253'  \
     --bulk-field='mri_swi_nii:20251'  \
     --bulk-field='mri_swi_dicom:20219'  \
     --bulk-field='mri_dmri_dicom:20218'  \
     --bulk-field='mri_dmri_nii:20250'  \
     --bulk-field='actigraphy_cwa:90001'  \
     --bulk-field='actigraphy_timeseries:90004'

# mri_t2_nii
screen -d -m bash -c "/ncf/sba01/simons_aging/slurm/download_bulk_id20253_download_all.bash"
# mri_rest_dicom
screen -d -m bash -c "/ncf/sba01/simons_aging/slurm/download_bulk_id20225_download_all.bash"
# mri_rest_nii
screen -d -m bash -c "/ncf/sba01/simons_aging/slurm/download_bulk_id20227_download_all.bash"
# actigraphy_cwa
screen -d -m bash -c "/ncf/sba01/simons_aging/slurm/download_bulk_id90001_download_all.bash"


cd /ncf/sba01/ukbAgingPipeline

base_dir=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline

for i in "rfmri_full_25 25750" "rfmri_full_100 25751" "rfmri_part_25 25752" "rfmri_part_100 25753" "rfmri_rsfa_25 25754" "rfmri_rsfa_100 25755"
do
    set -- $i
    python3 ${base_dir}/scripts/00b_download_bulk.py \
         --config=${base_dir}/config.json \
         --bulk-name=$1 \
         --bulk-id=$2 \
         --make-bulk-files \
         --slurm
done



for i in "rfmri_full_25 25750" "rfmri_full_100 25751"
do
    set -- $i
    python3 ${base_dir}/scripts/00b_download_bulk.py \
         --config=${base_dir}/config.json \
         --bulk-name=$1 \
         --bulk-id=$2 \
         --download-bulk-data \
         --slurm
done












