#!/bin/python

import os
import json
import subprocess
import argparse
import pandas as pd
import sys
import shutil

def main():
    example_text = '''example: 
     python test.py -t template/test.py'''

    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    parser.add_argument('--bulk-field', '-b', dest='bulk_field', action='append', nargs='+', required=True)
    parser.add_argument('--make-bulk-list', dest='make_bulk_list', action='store_true', default=False)
    parser.add_argument('--download-bulk-data', dest='download_bulk_data', action='store_true', default=False)
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    parser.add_argument('--slurm_partition', '-p', dest='slurm_partition', required=False, default=False)
    parser.add_argument('--singularity_container', dest='singularity_container', required=False, default=False)
    opt = parser.parse_args()

    print(opt.bulk_field)

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.config = '/ncf/sba01/ukbAgingPipeline/config.json'
        opt.bulk_field = [['rfmri_full_25:25750'], ['rfmri_full_100:25751'], ['rfmri_part_25:25752'], ['rfmri_part_100:25753'], ['rfmri_rsfa_25:25754'], ['rfmri_rsfa_100:25755']]
        opt.make_bulk_list = True
        opt.download_bulk_data = False
        opt.slurm = True
        opt.slurm_partition = 'short'
        opt.slurm_partition = 'ncf'
        opt.singularity_container = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/simons_ukb_aging_pipeline'
        opt.singularity_container = '/ncf/sba01/ukbAgingPipeline/simons_ukb_aging_pipeline'

    # parse args
    config_file     = opt.config
    bulk_field      = opt.bulk_field
    download_bulk_data = opt.download_bulk_data
    make_bulk_list  = opt.make_bulk_list
    slurm = opt.slurm
    slurm_partition = opt.slurm_partition
    singularity_container = opt.singularity_container

    if make_bulk_list == True and download_bulk_data == True:
        raise Exception("ERROR: '--download-bulk-data' and '--make-bulk-lists' both set to True. Only one may be set at a time. Run '--make-bulk-lists' first")

    # read config file
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # if submitting to slurm, import utilities
    if slurm == True:
        sys.path.append(os.path.join(config_json['repo_dir'], 'scripts'))
        sys.path.append(os.path.join(config_json['repo_dir'], 'scripts/utilities'))
        from utilities import submitSlurm, writeSlurm

    # folder to store the slurm commands
    slurm_dir = os.path.join(config_json['base_dir'], 'slurm')

    # point to the correct directories
    raw_dir  = os.path.join(config_json['base_dir'], 'data/ukb/raw')
    bulk_dir = os.path.join(config_json['base_dir'], 'data/ukb/bulk')

    # path to UKB encoding file
    enc_idx = 0
    ukb_enc  = config_json['ukb_encs'][enc_idx]['ukb_enc']
    enc_key  = config_json['ukb_encs'][enc_idx]['ukb_key']
    if ukb_enc[-4:] == '.enc':
        ukb_enc_file = ukb_enc.replace('.enc','.enc_ukb')
    elif ukb_enc[-8:] == '.enc_ukb':
        ukb_enc_file = ukb_enc
    else:
        raise Exception("CONFIG ERROR: 'ukb_enc' variable is mis-specified in the configuration file")

    if slurm == True and slurm_partition == False:
        raise Exception("ERROR: Must specify slurm partition name")

    # directory where all the operations will take place
    enc_dir  = '/'.join(ukb_enc.split('/')[:-1])
    enc_name = ukb_enc.split('/')[-1]
    enc_base = enc_name.split('.')[0]

    if make_bulk_list == True:
        print('Creating Bulk Files for Download')
        for download_field in bulk_field:
            bulk_name = download_field[0].split(':')[0]
            bulk_id   = download_field[0].split(':')[1]
            os.chdir(enc_dir)

            # where to write the data
            bulk_out_dir = os.path.join(bulk_dir, bulk_name)
            print('   Name:         {}'.format(bulk_name))
            print('   Field ID:     {}'.format(bulk_id))
            print('   Download Dir: {}'.format(bulk_out_dir))
            if not os.path.exists(bulk_out_dir):
                os.mkdir(bulk_out_dir)

            if slurm == True:
                print('....Submitting to Slurm....')

                # copy the ukbconv utility to the bulk download dir
                orig_ubconv = os.path.join(config_json['repo_dir'], 'external', 'ukbconv')
                shutil.copy(orig_ubconv, os.path.join(bulk_out_dir, 'ukbconv'))

                # symlink the ukb*enc_ukb file to bulk download dir
                src_enc_ukb  = os.path.join(raw_dir, '{}.enc_ukb'.format(enc_base))
                dest_enc_ukb = os.path.join(bulk_out_dir, '{}.enc_ukb'.format(enc_base))
                if not os.path.exists(dest_enc_ukb):
                    os.symlink(src_enc_ukb, dest_enc_ukb)

                # build command for slurm submission
                convert_cmd  = 'cd {}'.format(bulk_out_dir)
                download_cmd = ['./ukbconv', './{}.enc_ukb'.format(enc_base), 'bulk', '-s{}'.format(bulk_id), '-o./{}'.format(bulk_id)]
                submit_cmd   = '{}\n{}'.format(convert_cmd, ' '.join(download_cmd))

                # write slurm command file and submit to cluster
                slurm_file = os.path.join(slurm_dir, 'download_{}'.format(bulk_id))
                slurm_path = writeSlurm(slurm_file=slurm_file, partition=slurm_partition, cmd=submit_cmd, jobName=bulk_name, stime='6:00:00', n_gpu=None, nthreads=2, mem='10G')
                job_id     = submitSlurm(cmd=slurm_path)
                print('   Slurm Command: {}'.format(slurm_path))
                print('   Slurm Output:  {}'.format(slurm_path.replace('.txt','Out.txt')))
                print('   Slurm Job ID:  {}'.format(job_id))
                print('\n')

            else:
                write_file = os.path.join(slurm_dir, 'fetch_{}.txt'.format(bulk_id))
                with open(write_file, 'w') as f:
                    f.write(conv_cmd)
                print('Command to create bulk file written to here: \n{}'.format(write_file))

    if download_bulk_data == True:
        # if the bulk fields have already been created, then you can start the bulk download fetching
        print('Creating Bulk Files for Download')
        for download_field in bulk_field:
            bulk_name = download_field[0].split(':')[0]
            bulk_id   = download_field[0].split(':')[1]
            bulk_out_dir = os.path.join(bulk_dir, bulk_name)
            os.chdir(bulk_out_dir)

            # copy the ukbfetch utility to the bulk download dir
            orig_ukbconv = os.path.join(config_json['repo_dir'], 'external', 'ukbfetch')
            shutil.copy(orig_ukbconv, os.path.join(bulk_out_dir, 'ukbfetch'))
            key_name = enc_key.split('/')[-1]
            shutil.copy(enc_key, os.path.join(bulk_out_dir, key_name))

            # read file with subjects to pull
            bulk_file = os.path.join(bulk_out_dir, '{}.bulk'.format(bulk_id))
            bulk_df   = pd.read_csv(bulk_file, header=None, delim_whitespace=True)
            start = 1
            nrows = bulk_df.shape[0]
            job_array = []
            while start < nrows:
                print(start)
                slurm_file = os.path.join(slurm_dir, 'bulk_s{}_{}.bash'.format(start, bulk_id))
                job_array.append(slurm_file)

                fetch_cmd = f'''cd {bulk_out_dir}'''
                fetch_cmd = f'''{fetch_cmd}\n\n./ukbfetch \\\n {'-b./{}.bulk'.format(bulk_id)} \\\n -a./{key_name} \\\n -s{start} -m1000'''
                start = start + 1000

                # write command to slurm file
                if slurm == True:
                    slurm_fetch = os.path.join(slurm_dir, 'download_{}_{}'.format(start, bulk_id))
                    slurm_path = writeSlurm(slurm_fetch, slurm_partition, fetch_cmd, str(start), stime='6:00:00', n_gpu=None, nthreads=2, mem='8G')
                    print(slurm_path)
                    job_id = submitSlurm(slurm_path)

                else:
                    # write script to file
                    with open(slurm_file, 'w') as f:
                        f.write(fetch_cmd)
                    os.system('chmod 770 {}'.format(slurm_file))


if __name__ == "__main__":
    main()



