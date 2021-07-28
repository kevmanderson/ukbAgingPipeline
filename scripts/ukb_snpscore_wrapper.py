#!/bin/python

import os
import pandas as pd
import subprocess

# set up directories
base_dir  = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline'
out_dir   = os.path.join(base_dir, 'ref_files/gwas_sumstats/munged')
raw_dir   = os.path.join(base_dir, 'ref_files/gwas_sumstats/raw')
gwas_meta = pd.read_csv(os.path.join(base_dir, 'ref_files/simons_gwas_metadata_210723.csv'))
slurm_dir = '/gpfs/milgram/scratch/holmes/simons_aging/slurm'

# only do work on prepared gwas metadata
gwas_meta = gwas_meta.loc[gwas_meta['ready'] == True]

# ldsc script
ldsc_munge_py = os.path.join(base_dir, 'external/ldsc/munge_sumstats.py')


script_arr = []
for i,row in gwas_meta.iterrows():
    if pd.isna(row['munged_file']):
        print('skip')
        continue
    # i/o
    if pd.isna(row['unzipped_file']):
        sumstat_in  = os.path.join(raw_dir, '{}.txt'.format(row['munged_file']))
    else: 
        sumstat_in  = os.path.join(raw_dir, row['unzipped_file'])
    out_path    = os.path.join(out_dir, row['munged_file'])
    script_file = os.path.join(slurm_dir, 'munge_{}.txt'.format(row['munged_file']))


    # LDSC munge_stats command
    munge_base = f'''#!/bin/bash
source activate ldsc
cd {out_dir}
{ldsc_munge_py} \\
    --sumstats={sumstat_in} \\
    --out={ out_path } \\
    --a1={ row['A1'] } \\
    --a2={ row['A2'] } \\
    --p={ row['P'] } \\
    --signed-sumstats={ row['signed_sumstat'] } \\
    --snp={ row['SNP'] } \\'''


    # add total N number
    if pd.notna(row['N']):
        munge_base = f'''{munge_base}
--N={ row['N'] } \\'''

    elif pd.notna(row['Ncol']):
        munge_base = f'''{munge_base}
--N-col={ row['Ncol'] }'''

    # else add case/ctl N numbers
    elif pd.notna(row['Ncase']) and pd.notna(row['Nctl']):
        munge_base = f'''{munge_base}
--N-cas={ row['Ncase'] } \\
--N-con={ row['Nctl'] }'''

    # else define the case/ctl N column names
    elif pd.notna(row['Ncase_col']) and pd.notna(row['Nctl_col']):
        munge_base = f'''{munge_base}
--N-cas-col={ row['Ncase_col'] } \\
--N-con-col={ row['Nctl_col'] }'''
    
    # flag for ripke daner format
    if row['daner_n'] == True:
        munge_base = f'''{munge_base} --daner-n'''
    
    # write the LDSC munge command to disk
    with open(script_file, 'w') as f:
        f.writelines(munge_base)
    subprocess.call(['chmod', '755', script_file])
    script_arr.append(script_file)

# write text file of script paths to submite all at once through dsq/slurm
script_file = os.path.join(slurm_dir, 'munge_all_commands.txt')
with open(script_file, 'w') as f:
    f.writelines('\n'.join(script_arr))

dsq_cmd = f'''cd {slurm_dir}

dsq \\
--partition=psych_day \\
--mem-per-cpu=5gb \\
--cpus-per-task=2 \\
--nodes=1 \\
--job-file={script_file} \\
--output=/dev/null \\
--submit
'''
print(dsq_cmd)

    file_path = args.sumstat_in
    unzip_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/munged'


