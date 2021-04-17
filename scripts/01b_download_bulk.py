#!/bin/python

import os
import json
import subprocess
import argparse
import pandas as pd
import sys
import shutil
#sys.path.append('/scripts/utilities')

def submitSlurm(cmd, dependencies=None):
    if dependencies != None:
        # execute
        p = subprocess.Popen(['/usr/bin/sbatch', '--dependency=afterany:' + ':'.join(dependencies), cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    else:
        p = subprocess.Popen(['/usr/bin/sbatch', cmd], stdout=subprocess.PIPE)
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
        opt.bulk_field = [['rfmri_full_25:25750'], ['rfmri_full_100:25751'], ['rfmri_part_25:25752'], ['rfmri_part_100:25753'], ['rfmri_rsfa_25:25754'], ['rfmri_rsfa_100:25755']]
        opt.make_bulk_list = True
        opt.download_bulk_data = True
        opt.slurm = True
        opt.slurm_partition = 'short'
        opt.singularity_container = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/simons_ukb_aging_pipeline'

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

    if slurm == True:
        sys.path.append(os.path.join(config_json['repo_dir'], 'scripts'))
        sys.path.append(os.path.join(config_json['repo_dir'], 'scripts/utilities'))
        from utilities import submitSlurm, writeSlurm

    # create folder to store the slurm commands
    slurm_dir = os.path.join(config_json['base_dir'], 'slurm')
    if not os.path.exists(slurm_dir):
        os.mkdir(slurm_dir)

    # point to the correct directories
    raw_dir  = os.path.join(config_json['base_dir'], 'data/ukb/raw')
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

                cmd = 'singularity run {} \\\n'.format(singularity_container)
                cmd = f'''{cmd} python3 /scripts/\\\n'''.format(singularity_container)''

                singularity_container

                # copy the ukbconv utility to the bulk download dir
                shutil.copy('/app/ukbconv', os.path.join(bulk_out_dir, 'ukbconv'))

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
            # write script to file
            else:
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



