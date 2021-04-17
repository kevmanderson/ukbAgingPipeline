# UK Biobank Processing Pipeline

This repository processes UK Biobank data for integration within the Simons Brain Aging website.

### What This Repo Does:

1. Convert encoded UKB data into usable formats. 
2. Download and compiles "Bulk" UKB data.
3. Perform a modified PHESANT phenome-scan quantifying **AGE** effects. 
4. Compile UK Biobank variable metadata.
5. Dump all data into an SQL database. 


### Installation 

This code is meant to be run using Docker or Singularity. 

Feel free to clone this repository to modify/develop however you'd like. 

#### Singularity
```bash
singularity cache clean
singularity pull --name simons_ukb_aging_pipeline docker://kevinanderson/simons-ukb-aging-pipeline
```

#### Docker
```bash
docker pull kevinanderson/simons-bulk-rnaseq-pipeline
```


### Set-up Your Configuration File
```json
[
  {
    "repo_dir": "/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline",
    "base_dir": "/gpfs/milgram/project/holmes/kma52/buckner_aging",
    "ukb_enc": "/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc",
    "ukb_key": "/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.key"
  }
]
```
  
| Variable | Explanation |
| ------------- | ------------- |
| base_dir | Path to the primary project directory where results will be downloaded/written. The code will handle directory creation. This path is ideally empty, but its (probably) OK if not.  |
| repo_dir | Top-level path to the this code repository (i.e. where ```git clone``` put all the code). Should be different than ```base_dir```.  |
| ukb_enc | Full path to your raw + encoded UK Biobank data bucket, downloaded from the UKB AMS portal. Currently, the code only supports a single "omnibus" data bucket and doesn't combine fields that live in separate *enc_ukb* files. I suggest merging all your data into a single bucket using the UKB AMS tools...it will make your life easier in the long run anyways. | 
| ukb_key | Full path to your UKB key, required for decrypting the encoded ukb data| 


---

### Step 0: Prepare Directories
```bash
singularity run simons_ukb_aging_pipeline \
  python3 scripts/00_create_dirs.py \
    --config=./config.json
    
singularity run simons_ukb_aging_pipeline \
  python3 /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/00_create_dirs.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json
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

```bash
singularity run simons_ukb_aging_pipeline \
  python3 scripts/01_convert_ukbenc_data.py \
    --config=./config.json
    
singularity run simons_ukb_aging_pipeline \
  python3 /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/scripts/01_convert_ukbenc_data.py \
    --config=/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json
```

---

### Step 1 (Optional): Download Bulk MRI Data

Some of the neuroimaging phenotypes are not available in the ```*.enc_ukb``` file. These fields include resting-state "imaging derived phenotypes" (IDPs). We have to download and compile them separately. 

```bash
python3 00b_download_bulk.py \
         --config=FULL_DIR_PATH/config.json \
         --bulk-name="rfmri_full_25" \
         --bulk-id=25750 \
         --make-bulk-files
         # --slurm # use this option to submit download as SLURM job

# Let the above command finish executing before running this step. 
# second, actually download data
python3 00b_download_bulk.py \
         --config=FULL_DIR_PATH/config.json \
         --bulk-name="rfmri_full_25" \
         --bulk-id=25750 \
         --download-bulk-data
         # --slurm # use this option to submit download as SLURM job
```
| Bulk Name | Bulk Field ID | Description | 
| ------------- | ------------- | ------------- |
| rfmri_full_25 | 25750 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25750) |
| rfmri_full_100 | 25751 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25751) |
| rfmri_part_25 | 25752 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25752) |
| rfmri_part_100 | 25753 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25753) |
| rfmri_rsfa_25 | 25754 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25754) |
| rfmri_rsfa_100 | 25755 | [Jump to Showcase](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=25755) |

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



### Reference Files
See the README in **./ref_files/** directory for more detailed description of each file. 

* ```showcase.csv:``` Primary UKB metadata file with lots of phenotype information (see ref_files README)

* ```codings.csv:``` Data code mapping from the UK Biobank, matches numeric values to their textual meaning (e.g. 1=Yes).

* ```data-coding-ordinal-info.txt:``` Ordering/recoding information for ordinal variables.

* ```outcome-info.tsv:``` UKB metadata file from PHESANT. Has some unique fields not contained in showcase.csv. 




### Creating Synthetic Data




#### Contact
Developed by Kevin Anderson in Randy Buckner's Lab, in collaboration with Simon's Plasticity and Brain Aging


