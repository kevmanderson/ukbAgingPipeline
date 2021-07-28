#!/bin/python

import os
import glob
import subprocess
import numpy as np
import pandas as pd


# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging'
repo_dir = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline'

gwas_df = pd.read_csv(os.path.join(repo_dir, 'ref_files/simons_gwas_metadata_210723.csv'))


# convert bgen to bed file in prep for prsice
plink2   = '/gpfs/milgram/project/holmes/kma52/buckner_aging/external/plink2'
gene_dir = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/genetic/imputed_processed/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8'
for i in np.arange(1,23):

    bgen = '{}/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr{}.bgen'.format(gene_dir, i)
    sample = '{}/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr{}.sample'.format(gene_dir, i)
    out = '{}/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr{}'.format(gene_dir, i)
    script = '{}/bgen_to_bed_chr{}.bash'.format(gene_dir, i)
    out_file = '{}/bgen_to_bed_chr{}_output'.format(gene_dir, i)
    exclude = '{}/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr{}.rmdup.mismatch'.format(gene_dir, i)

    bgen_to_bed = f'''#!/bin/bash
    
{plink2} \\
--bgen {bgen} ref-first \\
--sample {sample} \\
--rm-dup \\
--make-bed \\
--exclude {exclude} \\
--out {out}
    '''
    print(bgen_to_bed)
    with open(script, 'w') as f:
        f.writelines(bgen_to_bed)
    subprocess.call(['chmod', '755', script])
    sbatch_cmd = 'sbatch --cpus-per-task=18 --mem-per-cpu=8GB --ntasks=1 --nodes=1 --partition=psych_day --time=360 --output={} {}'.format(out_file, script)
    os.system(sbatch_cmd)


# make a column with GWAS filename
gwas_df = gwas_df.loc[gwas_df['phenotype'].notna()]
gwas_df = gwas_df.loc[gwas_df['ready'] == True]
gwas_df['final_file'] = None
gwas_df.loc[gwas_df['unzipped_file'].notna(), 'final_file'] = gwas_df.loc[gwas_df['unzipped_file'].notna(), 'unzipped_file']
gwas_df.loc[gwas_df['unzipped_file'].isna(), 'final_file'] = gwas_df.loc[gwas_df['unzipped_file'].isna(), 'munged_file']


rscript      = '/gpfs/milgram/apps/hpc.rhel7/software/R/3.6.1-foss-2018b-X11-20180604/bin/Rscript'
prsice       = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/PRSice.R'
prsice_linux = os.path.join(repo_dir, 'external/PRSice_linux')
prs_dir      = os.path.join(base_dir, 'data/ukb/genetic/polygenic_scores')
sumstats_dir = os.path.join(repo_dir, 'ref_files/gwas_sumstats/raw')
gwas_data    = os.path.join(base_dir, 'data/ukb/genetic/imputed_processed/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr#')
extract      = os.path.join(base_dir, 'data/ukb/genetic/imputed_processed/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8/Alzheimers_2019_NatGen_Kunkle.valid')


# conduct PRSice for each GWAS
# -------
for i,gwas in gwas_df.iterrows():
    print(gwas['final_file'])

    # name of folder to read
    gwas_file = gwas['final_file']
    if '.txt' not in gwas_file:
        gwas_file = '{}.txt'.format(gwas_file)
    
    # folder for writing PRS score 
    folder_name  = gwas_file.split('.')[0]
    write_folder = os.path.join(prs_dir, folder_name)
    if not os.path.exists(write_folder):
        os.mkdir(write_folder)

    # gwas file
    sumstats = os.path.join(sumstats_dir, gwas_file)
    stat = gwas['signed_sumstat'].split(',')[0]
    null = gwas['signed_sumstat'].split(',')[1]

    prsice_cmd = f'''#!/bin/bash

{rscript} \\
{prsice} \\
--dir {write_folder} \\
--prsice {prsice_linux} \\
--base {sumstats} \\
--target {gwas_data} \\
--thread 4 \\
--A1 {gwas['A1']} --A2 {gwas['A2']} --stat {stat} \\
--snp {gwas['SNP']} --bp {gwas['POS']} --pvalue {gwas['P']} \\
--binary-target F \\
--fastscore \\
--bar-levels 1,0.05 \\
--no-regress \\
--extract {extract} \\
--out {write_folder}'''

    if null == '0':
        prsice_cmd = prsice_cmd + ' --beta'
    else: 
        prsice_cmd = prsice_cmd + ' --or'
    print(prsice_cmd)
    
    # write file
    out_file = os.path.join(prs_dir, 'run_prsice_{}_output.txt'.format(folder_name))
    cmd_file = os.path.join(prs_dir, 'run_prsice_{}.bash'.format(folder_name))
    with open(cmd_file, 'w') as f:
        f.writelines(prsice_cmd)
    subprocess.call(['chmod', '755', cmd_file])
    sbatch_cmd = 'sbatch --cpus-per-task=18 --mem-per-cpu=8GB --ntasks=1 --nodes=1 --partition=psych_day --time=360 --output={} {}'.format(out_file, cmd_file)
    os.system(sbatch_cmd)





# 
/gpfs/milgram/apps/hpc.rhel7/software/R/3.6.1-foss-2018b-X11-20180604/bin/Rscript \
    /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/PRSice.R \

/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/genetic/imputed_processed/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8


/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw/Alzheimers_2019_NatGen_Kunkle.txt

/gpfs/milgram/apps/hpc.rhel7/software/R/3.6.1-foss-2018b-X11-20180604/bin/Rscript \
    /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/PRSice.R \
    --dir /gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/genetic/polygenic_scores/Alzheimers_2019_NatGen_Kunkle \
    --prsice /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/PRSice_linux \
    --base /gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/ref_files/gwas_sumstats/raw/Alzheimers_2019_NatGen_Kunkle.txt \
    --target /gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/genetic/imputed_processed/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8/ukb_imp_maf0_01_hwe1e06_mind0_1_geno0_1_info0_8_chr# \
    --thread 4 \
    --stat BETA \
    --beta \
    --binary-target F \
    --fastscore --bar-levels 1,0.05 \
    --no-regress \
    --extract Alzheimers_2019_NatGen_Kunkle.valid \
    --out Alzheimers_2019_NatGen_Kunkle


