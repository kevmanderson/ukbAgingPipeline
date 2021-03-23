# UK Biobank Processing Pipeline



### Set-up Your Configuration File
```json
[
  {
    "base_dir": "/gpfs/milgram/project/holmes/kma52/buckner_aging",
    "ukb_enc": "/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc_ukb",
    "showcase_file": "/gpfs/milgram/project/holmes/kma52/buckner_aging/ref_files/showcase.csv",
    "codings_file": "/gpfs/milgram/project/holmes/kma52/buckner_aging/ref_files/codings.csv",
    "outcome_file": "/gpfs/milgram/project/holmes/kma52/buckner_aging/external/PHESANT/variable-info/outcome-info.tsv"
  }
]
```


### Step1: Prep UKB Data

We assume you've combined your data into a single omnibus data bucket, decrypted and stored wherever you plan to run your analyses (i.e. ```ukbXXXXX.enc_ukb```). 

Set the path to this file using the ```ukb_enc``` parameter in the config.json file. 

We convert the data using the ```r``` styled output for input to the web browser SQL database. We also generate a ```csv``` formatted version of the data to allow for input to a modified version of  [PHESANT](https://github.com/MRCIEU/PHESANT). 

```python
# Create formatted UKB data using this command (replace FULL_DIR_PATH with your path)
python 00_convert_ukbenc_data.py -c FULL_DIR_PATH/config.json
```


### Step 2: Prep UKB Metadata


### Reference Files
See the README in **./ref_files/** directory for more detailed description of each file. 

* ```showcase.csv:``` Primary UKB metadata file with lots of information (see ref_files README)

* ```codings.csv:``` Data code mapping from the UK Biobank, matches numeric values to their textual meaning (e.g. 1=Yes).

* ```data-coding-ordinal-info.txt:``` Ordering/recoding information for ordinal variables.

* ```outcome-info.tsv:``` UKB metadata file from PHESANT. Has some unique fields not contained in showcase.csv. 



### Creating Synthetic Data




#### Contact
Developed by Kevin Anderson in Randy Buckner's Lab, in collaboration with Simon's Plasticity and Brain Aging


