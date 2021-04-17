#!/bin/python

import os
import json
import subprocess
import argparse
import pandas as pd
import sys
sys.path.append('/scripts/utilities')
from utilities import submitSlurm, writeSlurm

def submitSlurm(cmd, dependencies=None):
    if dependencies != None:
        # execute
        p = subprocess.Popen(['sbatch', '--dependency=afterany:' + ':'.join(dependencies), cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    else:
        p = subprocess.Popen(['sbatch', cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    # get slurm job id to set up job submission dependency
    job_id = str(out).split(' ')[-1].replace("\\n'", '')
    return (job_id)


def writeSlurm(slurm_file, partition, cmd, jobName, stime='6:00:00', n_gpu=None, nthreads=None, mem=None):
    '''
    Submit batch job to Slurm-based cluster (SLURM)

    required arguments:
        - slurm_file        base filepath string for writing slurm command and output
                                (e.g. /gpfs/milgram/project/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/MRI/100024/slurm/100024_hello_world_cmd)
        - partition         short/long/scavenge
                                (e.g. short)
        - nthreads          up to 28 cpus on a nodemodel_performance
                                (e.g. 4)
        - cmd           full string of the command to be written/run
                            (e.g. module load Apps/FREESURFER/5.3.0\\n print('Hello World'))
    '''
    slurm_name = slurm_file + '_slurm.txt'
    slurm_out = slurm_file + '_slurmOut.txt'
    slurm_file = open(slurm_name, "w")
    slurm_file.write('#!/bin/bash\n')
    slurm_file.write('#SBATCH --partition=' + partition + '\n')
    slurm_file.write('#SBATCH --output=' + slurm_out + '\n')
    slurm_file.write('#SBATCH --nodes=1\n')
    if n_gpu != None:
        slurm_file.write('#SBATCH --gpus=' + str(n_gpu) + '\n')
    if mem != None:
        slurm_file.write('#SBATCH --mem=' + str(mem) + '\n')
    if nthreads != None:
        slurm_file.write('#SBATCH --ntasks=1 --cpus-per-task=' + str(nthreads) + '\n')
    slurm_file.write('#SBATCH --job-name=' + jobName + '\n')
    slurm_file.write('#SBATCH --time=' + stime + '\n')
    slurm_file.write(str(cmd))
    slurm_file.close()
    subprocess.call(['chmod', '0770', slurm_name])
    return (slurm_name)


def main():
    example_text = '''example: 
     python test.py -t template/test.py'''

    parser = argparse.ArgumentParser(epilog=example_text, add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', '-c', dest='config', required=True)
    parser.add_argument('--bulk-name', '-n', dest='bulk_name', required=True)
    parser.add_argument('--bulk-id', '-i', dest='bulk_id', required=True)
    parser.add_argument('--make-bulk-list', dest='make_bulk_list', action='store_true', default=False)
    parser.add_argument('--download-bulk-data', dest='download_bulk_data', action='store_true', default=False)
    parser.add_argument('--slurm', '-e', dest='slurm', action='store_true', default=False)
    opt = parser.parse_args()

    # read user configuration file
    tmp = False
    if tmp == True:
        parser = argparse.ArgumentParser()
        opt = parser.parse_args()
        opt.config = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
        opt.make_bulk_list = True
        opt.download_bulk_data = True
        opt.bulk_name = "rfmri_full_25"
        opt.bulk_id = 25750
        opt.slurm = True

    # parse args
    config_file = opt.config
    slurm = opt.slurm
    bulk_name = opt.bulk_name
    bulk_id = opt.bulk_id
    download_bulk_data = opt.download_bulk_data
    make_bulk_list = opt.make_bulk_list

    if make_bulk_list == True and download_bulk_data == True:
        raise Exception("User Error: '--download-bulk-data' and '--make-bulk-lists' both set to True. Only one may be set at a time. Run '--make-bulk-lists' first")

    # read config file
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]

    # create folder to store the slurm commands
    slurm_dir = os.path.join(config_json['base_dir'], 'slurm')
    if not os.path.exists(slurm_dir):
        os.mkdir(slurm_dir)

    bulk_dir = os.path.join(config_json['base_dir'], 'data/ukb/bulk')
    if not os.path.exists(bulk_dir):
        os.mkdir(bulk_dir)

    # path to UKB encoding file
    ukb_enc = config_json['ukb_enc']
    if ukb_enc[-4:] == '.enc':
        ukb_enc_file = ukb_enc.replace('.enc','.enc_ukb')
    elif ukb_enc[-8:] == '.enc_ukb':
        ukb_enc_file = ukb_enc
    else:
        raise Exception("Config Error: 'ukb_enc' variable is mis-specified in the configuration file")

    # directory where all the operations will take place
    enc_dir  = '/'.join(ukb_enc.split('/')[:-1])
    enc_name = ukb_enc.split('/')[-1]
    enc_base = enc_name.split('.')[0]

    # check to see if the ukbconv utility is downloaded
    ukb_conv_path  = os.path.join(enc_dir, 'ukbconv')
    ukb_fetch_path = os.path.join(enc_dir, 'ukbfetch')

    # dictionary of bulk fields to download
    bulk_dict = pd.Series({'rfmri_full_25': 25750,
                              'rfmri_full_100': 25751,
                              'rfmri_part_25': 25752,
                              'rfmri_part_100': 25753,
                              'rfmri_rsfa_25': 25754,
                              'rfmri_rsfa_100': 25755})

    if make_bulk_files == True:
        # We first have to create the bulk download files from the *enc.ukb file
        #for key in bulk_dict.keys():
        #field       = bulk_dict[key]

        conv_cmd = f'''cd {enc_dir}\n\n{ukb_conv_path} {enc_name} bulk -s{bulk_id} -o{bulk_id} -eencoding.ukb'''

        if slurm == True:
            slurm_fetch = os.path.join(slurm_dir, 'fetch_{}'.format(bulk_id))
            slurm_path  = writeSlurm(slurm_fetch, 'short', conv_cmd, bulk_id, stime='6:00:00', n_gpu=None, nthreads=2, mem='8G')
            job_id      = submitSlurm(slurm_path)
            print('Command to create bulk file written to here: \n{}'.format(slurm_path))
            print('Submitting job to slurm cluster')
        else:
            # write script to file
            write_file = os.path.join(slurm_dir, 'fetch_{}.txt'.format(bulk_id))
            with open(write_file, 'w') as f:
                f.write(conv_cmd)
            print('Command to create bulk file written to here: \n{}'.format(write_file))

    if download_bulk_data == True:
        # if the bulk fields have already been create, then you can start the bulk download fetching
        #for key in bulk_dict.keys():
            #field = bulk_dict[key]

        # create directory for all the output bulk files
        bulk_out_dir = os.path.join(enc_dir, 'bulk_{}'.format(bulk_id))
        if not os.path.exists(bulk_out_dir):
            os.mkdir(bulk_out_dir)

        # read file with subjects to pull
        bulk_file = os.path.join(enc_dir, '{}.bulk'.format(bulk_id))
        bulk_df = pd.read_csv(bulk_file, header=None, delim_whitespace=True)
        start = 1
        nrows = bulk_df.shape[0]
        job_array = []
        while start < nrows:
            print(start)
            slurm_file = os.path.join(slurm_dir, 'bulk_s{}_{}.bash'.format(start, bulk_id))
            job_array.append(slurm_file)

            fetch_cmd = f'''cd {bulk_out_dir}'''
            fetch_cmd = f'''{fetch_cmd}\n\n{ukb_fetch_path} \\\n {'-b../{}.bulk'.format(bulk_id)} \\\n -a../ukb40501.key \\\n -s{start} -m1000'''
            start = start + 1000

            # write command to slurm file
            if slurm == True:
                slurm_fetch = os.path.join(slurm_dir, 'download_{}_{}'.format(start, bulk_id))
                slurm_path = writeSlurm(slurm_fetch, 'short', fetch_cmd, str(start), stime='6:00:00', n_gpu=None, nthreads=2, mem='8G')
                job_id = submitSlurm(slurm_path)

            else:
                # write script to file
                with open(slurm_file, 'w') as f:
                    f.write(fetch_cmd)
                os.system('chmod 770 {}'.format(slurm_file))


if __name__ == "__main__":
    main()



