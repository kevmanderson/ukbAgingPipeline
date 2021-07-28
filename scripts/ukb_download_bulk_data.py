#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# Convert decrypted UKB data to readable formats

import os
import glob
import datetime
import calendar
import argparse
import logging
import numpy as np
import pandas as pd
import utilities.utilities as utilities

log = logging.getLogger('ukb')

def main():
    '''
    Before downloading bulk, data first need to make lists of subjecst that possess data

    Parameters
    ----------


    Returns
    -------

    '''
    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()

    parser = argparse.ArgumentParser()
    parser.add_argument('--ukb_fetch', dest='ukb_fetch', required=True, default=None)
    parser.add_argument('--enc_key', dest='enc_key', required=True, default=None)
    parser.add_argument('--bulk_field', dest='bulk_field', required=True, action='append', nargs='+', default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True)
    parser.add_argument('--slurm_dir', dest='slurm_dir', required=False)
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    args = parser.parse_args()
    print(args)

    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    args.bulk_field = [['rfmri_full_25:25750', 'rfmri_full_100:25751']]
    args.out_dir    = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/bulk'
    args.ukb_fetch  = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/external/ukb_tools/ukbfetch'
    args.slurm_dir  = '/gpfs/milgram/project/holmes/kma52/buckner_aging/slurm'
    args.enc_ukb    = '/gpfs/milgram/project/holmes/kma52/buckner_aging/data/ukb/raw/ukb40501.enc_ukb'

    log.info('RUNNING: MAKE_BULK_LIST')

    # iterate through all of the bulk fields to download
    for field in args.bulk_field:

        # get the descriptive bulk_name and ukb bulk_id
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]
        bulk_out_dir = os.path.join(args.out_dir, bulk_name)
        print(field)

        # read list of subjects to download
        bulk_file  = os.path.join(bulk_out_dir, '{}.bulk'.format(bulk_id))
        bulk_df    = pd.read_csv(bulk_file, header=None, delim_whitespace=True)

        # check which files have already been downloaded
        downloaded_paths  = glob.glob(os.path.join(bulk_out_dir, '*{}*'.format(bulk_id)))
        downloaded_files  = set([x.split('/')[-1] for x in downloaded_paths])

        # if no files downloaded
        if len(downloaded_files) > 0:
            append = '{}'.format(list(downloaded_files)[0].split('.')[1])
        else:
            append = 'txt'

        # create list of filenames
        total_files = ['{}_{}.{}'.format(i[0], i[1], append) for i in zip(bulk_df[0], bulk_df[1])]

        # find difference between downloaded and to-be-downloaded
        files_to_download = list(set(total_files) - set(downloaded_files))

        if len(files_to_download) > 0:

            download_bulk_file = bulk_file.replace('.bulk', '_use.bulk')
            subs_to_download   = list(set([x.split('_')[0] for x in files_to_download]))
            bulk_download_df   = bulk_df.loc[bulk_df[0].astype(str).isin(subs_to_download)]
            bulk_download_df.to_csv(download_bulk_file, sep='\t', header=None, index=None)

            fetch_cmd = f'''cd {bulk_out_dir}'''
            slurm_file = os.path.join(write_dir, 'download_bulk_id{}_download_all.bash'.format(bulk_id))

            chunk_size = 1000
            start = 1
            nrows = bulk_download_df.shape[0]
            job_array = []
            while start < nrows:
                #slurm_file = os.path.join(write_dir, 'download_bulk_id{}_s{}.bash'.format(bulk_id, start))
                job_array.append(slurm_file)

                # fetch command
                fetch_cmd = f'''{args.ukb_fetch}\n./ukbfetch {'-b./{}_use.bulk'.format(bulk_id)} -a./{key_file} -s{start} -m{chunk_size}'''
                start = start + chunk_size

            if args.slurm == True:
                slurm_path = writeSlurm(slurm_file,
                                        slurm_partition,
                                        fetch_cmd,
                                        str(start),
                                        stime='6:00:00',
                                        n_gpu=None,
                                        nthreads=2,
                                        mem='8G')
                job_id = submitSlurm(slurm_path)
            else:
                utilities.write_cmd(cmd=fetch_cmd, write_file=slurm_file)




        enc_name = args.enc_ukb.split('/')[-1]

        # build command for slurm submission
        convert_cmd  = 'cd {}'.format(bulk_out_dir)
        download_cmd = [args.ukb_fetch, './{}'.format(enc_name), 'bulk', '-s{}'.format(bulk_id), '-o./{}'.format(bulk_id)]
        submit_cmd   = '{}\n{}'.format(convert_cmd, ' '.join(download_cmd))
        print(submit_cmd)

        write_file = os.path.join(args.slurm_dir, 'make_bulk_sublist_{}.bash'.format(bulk_id))
        utilities.write_cmd(cmd=submit_cmd, write_file=write_file)


if __name__ == "__main__":
    main()




