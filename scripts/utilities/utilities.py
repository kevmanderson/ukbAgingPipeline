#!/bin/python

import os
import subprocess
import sys



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
    Submit batch job to Yale Milgram cluster (SLURM)

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


def submitSlurmArrays(job_array, jobs_in_batch, slurm_base, ncpus, partition, jobName, stime='6:00:00', mem=None):
    batch_list = list(chunks(job_array, jobs_in_batch))
    # write chunked batch list
    for batch_num,batch_files in enumerate(batch_list):
        print(batch_num)
        slurm_batch = os.path.join(slurm_base + '_batch_{}.txt'.format(batch_num))
        slurm_file  = open(slurm_batch, "w")
        slurm_file.write('#!/bin/bash\n')
        for write_file in batch_files:
            slurm_file.write(write_file + '\n')
        slurm_file.close()
        subprocess.call(['chmod', '0770', slurm_batch])

        slurm_dir  = '/'.join(slurm_base.split('/')[:-1])
        batch_file = os.path.join(slurm_base + '_batch_{}.sh'.format(batch_num))
        dsq_cmd = f'''/gpfs/milgram/apps/hpc.rhel7/software/dSQ/1.05/dsq --partition={partition}'''
        dsq_cmd = f'''{dsq_cmd} --cpus-per-task={ncpus}'''
        dsq_cmd = f'''{dsq_cmd} -J {jobName}'''
        dsq_cmd = f'''{dsq_cmd} --batch-file {batch_file}'''
        dsq_cmd = f'''{dsq_cmd} --job-file {slurm_batch}'''
        dsq_cmd = f'''{dsq_cmd} --status-dir {slurm_dir}'''
        dsq_cmd = f'''{dsq_cmd} --submit -t {stime}'''
        os.system(dsq_cmd)
















