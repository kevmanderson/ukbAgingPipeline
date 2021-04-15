#!/bin/bash

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












