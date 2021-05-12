#!/bin/python

import os
import subprocess
import sys
import itertools



def read_funpack_categories(cat_file):
    '''Read a funpack formatted category.tsv file'''
    with open(cat_file, 'r') as f:
        hdr = f.readline()
        line = f.readlines()
    field_string = [l.split('\t')[-1].replace('\n', '').replace('"', '') for l in line]
    field_string = reduce(lambda x, y: x + ',' + y, field_string)
    field_list = field_string.split(',')

    parsed_fields_l = [parseMatlabRange(l) for l in field_list]
    parsed_fields = reduce(lambda x, y: x + y, parsed_fields_l)

    return parsed_fields


def write_cmd(cmd, write_file, print_path=False):
    '''Write single or multi-line command to file'''
    with open(write_file, 'w') as f:
        f.write(cmd)
    os.chmod(write_file, 755)
    if print_path == True:
        'Writing command: {}'.format(print_path)

def make_symlink(src_path, dest_path):
    '''Symlink file. Remove and recreate if it already exists'''
    if src_path != dest_path:
        if os.path.exists(dest_path):
            print('Removing previously created symlinked file: {}'.format(dest_path))
            os.remove(dest_path)
        print('Symlink: {}'.format(dest_path))
        os.symlink(src_path, dest_path)
    else:
        print('src file is the same as dest file, continuing.')


def make_dir(create_path):
    '''Create a path if it doesnt exist, but also print feedback'''
    if not os.path.exists(create_path):
        os.mkdir(create_path)
        print('Directory created: {}'.format(create_path))
    else:
        print('Directory exists: {}'.format(create_path))


def create_directories(root_dir):
    '''
    Create project directories

    :param root_dir: directory where data will be written
    :return:
    '''

    make_dir(root_dir)
    dir_list = ['data', 'slurm', 'data/ukb',  'data/ukb/external', 'data/ukb/bulk', 'data/ukb/raw',
                      'data/ukb/genetic', 'data/ukb/genetic/genotyped', 'data/ukb/genetic/imputed']
    for create_path in dir_list:
        make_dir(os.path.join(root_dir, create_path))



def submit_slurm(cmd, dependencies=None):
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


def write_slurm(slurm_file, partition, cmd, jobName, stime='6:00:00', n_gpu=None, nthreads=None, mem=None):
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



def parseMatlabRange(r):
    """Parses a string containing a MATLAB-style ``start:stop`` or
    ``start:step:stop`` range, where the ``stop`` is inclusive).

    :arg r:   String containing MATLAB_style range.
    :returns: List of integers in the fully expanded range.
    """
    elems = [int(e) for e in r.split(':')]

    if len(elems) == 3:
        start, step, stop = elems
        if   step > 0: stop += 1
        elif step < 0: stop -= 1

    elif len(elems) == 2:
        start, stop  = elems
        stop        += 1
        step         = 1
    elif len(elems) == 1:
        start = elems[0]
        stop  = start + 1
        step  = 1
    else:
        raise ValueError('Invalid range string: {}'.format(r))

    return list(range(start, stop, step))













