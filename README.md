# UK Biobank Processing Pipeline

Author: Kevin Anderson (kevinanderson@fas.harvard.edu), Post-Doc, Harvard University.

This repository was created in the Buckner Laboratory.

It documents how we process UK Biobank data for use in a related interactive browser [insert link]

Running this code will require a bit of unavoidable babysitting, since most stages need to be run sequentially. For instance, bulk MRI data has to first be downloaded from the UKB before it can be combined and processed. 


### Table of Contents

| | Process | Link | Explanation |
| ------------- | ------------- | ------------- | ------------- |
| 1 | Installation | [details](#installation) |  How to install this repository | 
| 2 | Configuration | [details](#required-configuration-file) | Configuration file with required data paths | 
| 3 | Directory Structure | [details](#directory-structure) | Create project directories that will be populated with data  | 
| 4 | Decrypt UKB Data | [details](#decrypt-ukb-data) | Decrypt `ukb*.enc` to `ukb*.enc_ukb` format | 
| 5 |Convert UKB Data | [details](#decrypt-ukb-data) | Converts `ukb*.enc_ukb` files to usable csvs and tables | 
| 6 | Download Bulk MRI Data | [details](#download-bulk-data) | Bulk MRI data download (conducted in 3 sequential steps) | 
| 7 | PHESANT Pipeline | [details](#phesant-pipeline) | Phenome-wide regression against AGE. PHESANT also minimally preprocesses UKB phenotype data | 
| 8 | Genetic Pipeline | [details](#genetic-pipeline) | Download and preprocess UKB SNP data | 




### Installation 

Feel free to clone this repository to modify and develop however you'd like. 

#### Option 1: Clone repository
```bash
# use this option if you'd like to modify/develop any code
git clone https://github.com/kevmanderson/buckner-lab-ukb-pipeline.git

# build local container from within the code repository
docker -t buckner_lab_ukb_pipeline build .
```

#### Option 2: Singularity
Singularity is usually preferred in cluster environments. 
```bash
# optional, clear prior singularity images
singularity cache clean

# Pull the container from DockerHub (w/ Singularity)
singularity pull --name buckner_lab_ukb_pipeline docker://kevinanderson/buckner-lab-ukb-pipeline
singularity shell buckner_lab_ukb_pipeline
```

#### Option 3: Docker
```bash
# Pull the container from DockerHub (w/ Docker)
docker pull kevinanderson/buckner-lab-ukb-pipeline
```


## Required Configuration File

This pipeline requires you to specify some data paths and directories specific to your environment. 

See below for an example: 

```json
[
  {
    "repo_dir": "/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline",
    "base_dir": "/gpfs/milgram/project/holmes/kma52/buckner_aging",
    "ukb_encs" : [
      {
        "ukb_enc": "/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc",
        "ukb_key": "/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.key"
      }],
    "genotyped_data": "/gpfs/milgram/data/UKB/REPOSITORY/GWAS/NON_IMPUTED/ukb_cal_chr*_v2",
    "imputed_data": [
      {
        "bgen": "/gpfs/milgram/data/UKB/ukb_snp/ukb_imp_chr*_v3",
        "sample": "/gpfs/milgram/data/UKB/ukb_snp/ukb25163_imp_chr*_v3_s487324"
      }]
  }
]
```
  
| Variable | Explanation |
| ------------- | ------------- |
| base_dir | Path to the primary project directory where results will be downloaded/written. The code will handle directory creation. This path is empty, but its (probably) OK if not. Path should have a healthy (>100GB) amount of disk storage. |
| repo_dir | Top-level path to the this code repository (i.e. where ```git clone``` put all the code). Should be different than ```base_dir```.  |
| ukb_encs | List of ukb_enc files and their keys. UKB data are often split into multiple "buckets". This pipeline with read and aggregate separate data buckets into a single source. |
| ukb_enc | Full path to your raw UK Biobank data bucket, downloaded from the UKB AMS portal. I suggest merging all your data into a single bucket using the UKB AMS tools...it will make your life easier in the long run anyways. | 
| ukb_key | Full path to your UKB key, required for decrypting the encoded ukb data |
| genotyped_data | [OPTIONAL] Full path to already downloaded UKB genotype (not imputed bgen) files. Replace the chromosome number with a wildcard "*" |
| imputed_data | [OPTIONAL] Full path to already downloaded UKB imputed bgen + sample files (not genotyped plink). |
| bgen | [OPTIONAL] Path to bgen files (replace chromosome number with "*") |
| sample | [OPTIONAL] Path to sample files (replace chromosome number with "*") |

---
---
## Directory Structure

Create the directory structure for downloading and processing data.

```
. "base_dir"
└─── external
└─── slurm
└─── data
     └─── ukb
          └─── bulk
          └─── raw
          └─── output
              └─── sqlite_files
          └─── genotypes
              └─── genotyped
              └─── imputed
```
---
---

## Decrypt UKB Data

The first step is to decrypt the *enc files provided by the UK Biobank using the ```ukbunpack``` utility

See the [UKB Documentation](https://biobank.ctsu.ox.ac.uk/~bbdatan/Accessing_UKB_data_v2.3.pdf) for details. 

```bash
# Example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
  --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
  --stage='decrypt'
```

**input/output**
```bash
# Example Input: 
${repo_dir}/data/ukb/raw/ukb40501.enc

# Example Output: 
${repo_dir}/data/ukb/raw/ukb40501.enc_ukb
```
---
---

## Convert UKB Data

Once the data buckets have been converted to ```ukb*.enc_ukb```, we will convert from binary to human readable csv files using the ```ukbconv``` command.

This code separately creates a smaller csv file containing covariates required for PHESANT regression. These UKB fields are:
* 31 - sex
* 54 - assessment center
* 21003 - age at assessment center
* 22009 - genetic principle components


```bash
# Example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
  --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
  --stage='convert'
```  

`input/output`
```bash
# Example Input: 
${repo_dir}/data/ukb/raw/ukb40501.enc_ukb

# Example Output: 
${repo_dir}/data/ukb/raw/ukb40501.r
${repo_dir}/data/ukb/raw/ukb40501.txt
${repo_dir}/data/ukb/raw/ukb40501.csv
```
---
---

## Download Bulk Data

These stage will download bulk fields. We document downloads for RSFC/RSFA data, but this works equally well for any bulk UKB field. 

Acquiring the data requires 3 sequantial steps:  
(1) PREP:     Make download lists of available data  
(2) DOWNLOAD: Download individual bulk fields  
(3) FORMAT:   Compile individual files into single csv dataframe

**01: PREP**
```bash
# Example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='make_bulk_list' \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755'  
```

#### output
```bash
# Example Output: 
${repo_dir}/data/ukb/bulk/rfmri_full_25/25750.bulk
```


**02: DOWNLOAD**
```bash
# Example command

cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='download_bulk' \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755'  
```

#### input/output
```bash
# Example Input:
${repo_dir}/data/ukb/bulk/rfmri_full_25/25750.bulk

# Example Output: 
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id1}_25750_2_0.txt
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id2}_25750_2_0.txt
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id3}_25750_2_0.txt
```

**03: FORMAT**
```bash
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='bulk_to_csv' \
    --bulk-field='rfmri_full_25:25750'  \
    --bulk-field='rfmri_full_100:25751'  \
    --bulk-field='rfmri_part_25:25752'  \
    --bulk-field='rfmri_part_100:25753'  \
    --bulk-field='rfmri_rsfa_25:25754'  \
    --bulk-field='rfmri_rsfa_100:25755'  
```
  
#### input/output
```bash
# Example Input: 
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id1}_25750_2_0.txt
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id2}_25750_2_0.txt
${repo_dir}/data/ukb/bulk/rfmri_full_25/${ukb_id3}_25750_2_0.txt

# Example Output:
${repo_dir}/data/ukb/raw/rfmri_full_25/bulk_25750_2_0.csv
```
---
---

## PHESANT Pipeline

Data are next run through a custom version of PHESANT, once data are downloaded and formatted into csv files from above.

The PHESANT pipeline runs an age-regression for each phenotype, and saves the preprocessed data. 

Sequential steps:  
(1) PREPARE: Put the csv files in the correct format for PHESANT  
(2) RUN: Run PHESANT regression and save preprocessed data  

__01: PREPARE__
```bash
# example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

# Run this for each data bucket 
python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='prep_data_for_phesant' \
    --phesant-covar-csv-list='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_covars.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25750_3.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25751_3.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25752_3.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25753_3.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25754_3.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_2.csv' \
    --phesant-csv='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/bulk_25755_3.csv'


# temporary code for second data bucket    
python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='prep_data_for_phesant' \
    --phesant-covar-csv-list='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_covars.csv' \
    --phesant-data-csv-list='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb43410.csv' 

```

If you are running this with bulk RSFA & RSFC data

```bash
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

# iterate over all the bulk mri filenames
array=( bulk_25750_2 bulk_25750_3 bulk_25751_2 bulk_25751_3 bulk_25752_2 bulk_25752_3 bulk_25753_2 bulk_25753_3 bulk_25754_2 bulk_25754_3 bulk_25755_2 bulk_25755_3 ) 
for var in "${array[@]}"
do
python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='prep_data_for_phesant' \
    --phesant-covar-csv-list='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_covars.csv' \
    --phesant-data-csv-list='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/'${var}'.csv' 
done
```


#### run PHESANT

__01: RUN__
```bash
# example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='run_phesant' \
    --phesant-visits='0;1;2;3' \
    --phesant-phenofile='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit*_regress.csv' \
    --phesant-phenofile='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501_phesant_visit*_process.csv' \
    --phesant-phenofile='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb43410_phesant_visit*_regress.csv' \
    --phesant-phenofile='/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb43410_phesant_visit*_process.csv' \
    --slurm \
    --slurm_partition='short'
```  
    


## Genetic Pipeline

---
---
Genetic preprocessing is conducted on imputed UKB genetic data. 

```bash
# example command
cd /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline
source ukb_venv/bin/activate

python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='download_genetic' 
    
python ./scripts/main.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
    --stage='snp_preprocess' \
    --maf=0.01 \
    --hwe=1e-6 \
    --mind=0.1 \
    --geno=0.1 \
    --info=0.8 \
    --slurm
    
```  




# Example command (obviously, define $repo_dir variable yourself)
cd ${repo_dir}
singularity run simons_ukb_aging_pipeline \
  python3 /scripts/00_create_dirs.py \
    --config=./config.json

# if modifying code, use full path to the script
cd ${repo_dir}
singularity run simons_ukb_aging_pipeline \
  python3 ${repo_dir}/scripts/00_create_dirs.py \
    --config=${repo_dir}/config.json
   

# TMP -- harvard 
singularity run simons_ukb_aging_pipeline \
  python3 /ncf/sba01/ukbAgingPipeline/scripts/00_create_dirs.py \
    --config=/ncf/sba01/ukbAgingPipeline/harvard_config.json
    
# TMP -- yale    
singularity run simons_ukb_aging_pipeline \
  python3 /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/00_create_dirs.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json
```
---

### Step 1: Prep UKB Raw Data

Next, decrypt the ```ukb_enc``` file and convert to the proper formats.

We assume you've merged your data into a single omnibus data bucket (i.e. ```ukb*****.enc```). 

Set the path to this file using the ```ukb_enc``` parameter in the ```config.json``` file.

We decrypt the ```*.enc``` file using the [ukbunpack](https://biobank.ndph.ox.ac.uk/showcase/download.cgi) utility. 

We convert the ```*.enc_ukb``` file to readable formats using the [ukbconv](https://biobank.ndph.ox.ac.uk/showcase/download.cgi) utility. 

Calling the script below will result in two outputs:   
    1. ```r```: Primary format used for creating the SQL database.  
    2. ```csv```: Required for running a modified version of [PHESANT](https://github.com/MRCIEU/PHESANT). 

Execute the data preparation step with the following command:

*N.B.* This will take a __30-90__ minutes...time for coffee.

```bash
singularity run simons_ukb_aging_pipeline \
  python3 scripts/01_convert_ukbenc_data.py \
    --config=./config.json

# harvard 
singularity run simons_ukb_aging_pipeline \
  python3 /ncf/sba01/ukbAgingPipeline/scripts/01_convert_ukbenc_data.py \
    --config=/ncf/sba01/ukbAgingPipeline/config.json
    
# yale    
singularity run simons_ukb_aging_pipeline \
  python3 /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/01_convert_ukbenc_data.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json
```

---

### Step 1b (Optional): Download Bulk MRI Data

Some of the neuroimaging phenotypes are not available in the ```*.enc_ukb``` file. These fields include resting-state "imaging derived phenotypes" (IDPs). We have to download and compile them separately. 

A) First, create *.bulk files listing the bulk data to download. 
```bash
# harvard
cd /ncf/sba01/ukbAgingPipeline
python3 ./scripts/01b_download_bulk.py \
         --config=/ncf/sba01/ukbAgingPipeline/config.json \
         --bulk-field='mri_t1_nii:20252'  \
         --make-bulk-list \
         --slurm \
         --slurm_partition='ncf'
         
cd /ncf/sba01/ukbAgingPipeline
python3 ./scripts/01b_download_bulk.py \
         --config=/ncf/sba01/ukbAgingPipeline/config.json \
         --bulk-field='rfmri_full_25:25750'  \
         --bulk-field='rfmri_full_100:25751'  \
         --bulk-field='rfmri_part_25:25752'  \
         --bulk-field='rfmri_part_100:25753'  \
         --bulk-field='rfmri_rsfa_25:25754'  \
         --bulk-field='rfmri_rsfa_100:25755'  \
         --bulk-field='mri_rest_dicom:20225'  \
         --bulk-field='mri_rest_nii:20227'  \
         --bulk-field='mri_t1_nii:20252'  \
         --bulk-field='mri_t2_nii:20253'  \
         --bulk-field='mri_swi_nii:20251'  \
         --bulk-field='mri_swi_dicom:20219'  \
         --bulk-field='mri_dmri_dicom:20218'  \
         --bulk-field='mri_dmri_nii:20250'  \
         --bulk-field='actigraphy_cwa:90001'  \
         --bulk-field='actigraphy_timeseries:90004' \
         --make-bulk-list \
         --slurm \
         --slurm_partition='ncf'
```

B) Second, actually download the bulk data. 
```bash
# harvard
python3 ./scripts/01b_download_bulk.py \
         --config=/ncf/sba01/ukbAgingPipeline/config.json \
         --bulk-field='actigraphy_cwa:90001'  \
         --download-bulk-data \
         --slurm \
         --slurm_partition='ncf'

# yale
python3 ./scripts/01b_download_bulk.py \
         --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
         --bulk-field='rfmri_full_25:25750'  \
         --bulk-field='rfmri_full_100:25751'  \
         --bulk-field='rfmri_part_25:25752'  \
         --bulk-field='rfmri_part_100:25753'  \
         --bulk-field='rfmri_rsfa_25:25754'  \
         --bulk-field='rfmri_rsfa_100:25755'  \
         --download-bulk-data \
         --slurm \
         --slurm_partition='short'
```

C) Once bulk MRI data have been downloaded, read and compile them into dataframes.

```bash
singularity run simons_ukb_aging_pipeline \
python3 ./scripts/01c_compile_bulk_data.py \
         --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/yale_config.json \
         --bulk-field='rfmri_full_25:25750'  \
         --bulk-field='rfmri_full_100:25751'  \
         --bulk-field='rfmri_part_25:25752'  \
         --bulk-field='rfmri_part_100:25753'  \
         --bulk-field='rfmri_rsfa_25:25754'  \
         --bulk-field='rfmri_rsfa_100:25755'
```

#### Neuroimaging Bulk Fields

#### TODO: finish table/field descriptions

| Bulk Name | Bulk Field ID | Description | 
| ------------- | ------------- | ------------- |
| rfmri_full_25 | 25750 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25750) |
| rfmri_full_100 | 25751 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25751) |
| rfmri_part_25 | 25752 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25752) |
| rfmri_part_100 | 25753 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25753) |
| rfmri_rsfa_25 | 25754 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25754) |
| rfmri_rsfa_100 | 25755 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25755) |

---

### Step 2: Download Genetic Data

```bash

python3 ./scripts/02_download_genetic.py \
         --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json \
         --bulk-field='rfmri_full_25:25750'  \
         --bulk-field='rfmri_full_100:25751'  \
         --bulk-field='rfmri_part_25:25752'  \
         --bulk-field='rfmri_part_100:25753'  \
         --bulk-field='rfmri_rsfa_25:25754'  \
         --bulk-field='rfmri_rsfa_100:25755'

```

---

### Step 3: Prep UKB Metadata

The next step is to combine UKB variable information into a single metadata file (```ukbMetaData.csv```). 

We join data from the [UKB showcase](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide) with some curated fields from [PHESANT](https://github.com/MRCIEU/PHESANT/tree/master/variable-info). 

The final output of this step is the ```ukbMetaData.csv``` dataframe:

| Column Name | Type | Description |  
| ------------- | ------------- | ------------- |
| FieldID  | Int | UK Biobank Field. Each variable has a unique numeric id that exactly matches the  [UKB Showcase](https://biobank.ndph.ox.ac.uk/showcase/browse.cgi) |
| Category  | Int | Numeric ID for category. Each category contains some number of related UKB fields (e.g. 152=Process Durations, 100027=Fluid Int items) |
| Path | String | UKB fields are organized into a hierarchy. This variables contains the full categorical ontology (e.g. "UK Biobank Assessment Centre > Procedural metrics > Process durations") |
| Field | String | Text description of the UKB field. |
| ValueType | String | 1 of 5 main data types in the UKB dataset (i.e. Integer, Categorical single, Date, Continuous, Categorical multiple)  |
| DataCoding | Int | Many UKB fields require recoding from integer to string. This provides the unique ID for variable recoding information, which can be referenced to the ```datacodes``` SQL table. |
| Excluded | String | PHESANT field. Hand marks variables that should be excluded from phenotypic regressions (e.g. date fields). See [PHESANT](https://github.com/MRCIEU/PHESANT) |
| CatMultIndicatorFields | String | See [PHESANT](https://github.com/MRCIEU/PHESANT). |
| CatSingleToCatMult | String | See [PHESANT](https://github.com/MRCIEU/PHESANT). |
| DateConvert | String | This is all NAs for now... |
| Participants | Int | Number of subjects with at least one instance of this variable. |
| Items | Int | Number of unique instances of this variable (across all visits). |
| Units | String | Variable unit of measurement (e.g. seconds, mm3). |
| Strata | String | Auxiliary, Primary, Derived, Supporting |
| Instances | Int | Number of times this variable was measured. |
| Array | Int | TODO |
| Notes | String | Text based miscellaneous notes about the variable |
| Link | String | URL link to the variable info page | 
| phenoCategory | String | Most specific category for this variable. |
| covarType | String | Which initial covariates to use in the regression. |  


This step will also create the ```covariateTable.csv``` file used for preselecting confounders used in the regression. 

---

### Step 4: Create UKB SQLite file

---

### Step 5: Genetic Preprocessing

Genetic preprocessing is conducted on imputed UKB genetic data. 

#####Stage 1: Initial Filter
A. Variant filtering to retrain.  
2. Bi-allelic variants.  
3. SNPs with imputed accuracies (i.e. INFO) > 0.60.  

B. Individual filtering to retain individuals with:  
1. European Ancestry.
2. No Sex Aneuploidy.
3. Genetic sex matching reported sex. 
4. Inclusion in UKB Kinship estimation. 

##### Stage 2: Combined Filtered Data

bgen files are originally split by chromsome, we combined them into a single file to allow for QC with plink. 

#### Stage 3: Plink QC


```bash
# Option 1: run the command locally on each chromosome
for chrom in {1..22};
do
    echo $chrom
    python3 ./scripts/05_ukb_genetic_pipe.py \
             --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json \
             --chrom=${chrom} \
             --info=0.06 \
             --filter_bgen
done
    
    
# Option 2: Run each chromosome filtering in parallel using slurm
slurm_dir=/gpfs/milgram/project/holmes/kma52/buckner_aging/slurm
for chrom in {1..22};
do
    slurm_file=${slurm_dir}/snp_qc_filter_bgen_chr${chrom}.txt
    slurm_out_file=${slurm_dir}/snp_qc_filter_bgen_chr${chrom}_out.txt
    echo $slurm_file
    # write singularity command to slurm submission file
cat << EOF > ${slurm_file}
#!/bin/bash
#SBATCH --partition=short
#SBATCH --output=${slurm_out_file}
#SBATCH --nodes=1
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --job-name=$chrom
#SBATCH --time=06:00:00

singularity run simons_ukb_aging_pipeline \\
python3 /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/05_ukb_genetic_pipe.py \\
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json \\
    --chrom=${chrom} \\
    --maf=0.01 \\
    --hwe=1e-6 \\
    --mind=0.02 \\
    --geno=0.02 \\
    --info=0.6 \\
    --filter_bgen
echo 'done'
EOF
    # submit job to cluster
    sbatch < ${slurm_file}
done
    
    
```

---

### Reference Files
See the README in **./ref_files/** directory for more detailed description of each file. 

* ```showcase.csv:``` Primary UKB metadata file with lots of phenotype information (see ref_files README)

* ```codings.csv:``` Data code mapping from the UK Biobank, matches numeric values to their textual meaning (e.g. 1=Yes).

* ```data-coding-ordinal-info.txt:``` Ordering/recoding information for ordinal variables.

* ```outcome-info.tsv:``` UKB metadata file from PHESANT. Has some unique fields not contained in showcase.csv. 




### Creating Synthetic Data




#### Contact
Developed by Kevin Anderson in Randy Buckner's Lab, in collaboration with Simon's Plasticity and Brain Aging


