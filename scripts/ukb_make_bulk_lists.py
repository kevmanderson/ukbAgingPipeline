#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
# Convert decrypted UKB data to readable formats

import os
import datetime
import calendar
import argparse
import logging
import numpy as np
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
    parser.add_argument('--ukb_conv', dest='ukb_conv', required=True, default=None)
    parser.add_argument('--enc_ukb', dest='enc_ukb', required=True, default=None)
    parser.add_argument('--enc_key', dest='enc_ukb', required=True, default=None)
    parser.add_argument('--bulk_field', dest='bulk_field', required=True, action='append', nargs='+', default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True)
    parser.add_argument('--slurm_dir', dest='slurm_dir', required=False)
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    args = parser.parse_args()
    print(args)

    log.info('RUNNING: MAKE_BULK_LIST')

    # iterate through all of the bulk fields to download
    for field in args.bulk_field:
        # get the descriptive bulk_name and ukb bulk_id
        bulk_name = field[0].split(':')[0]
        bulk_id   = field[0].split(':')[1]
        print(field)

        # where the bulk data will be downloaded
        bulk_out_dir = os.path.join(args.out_dir, bulk_name)
        utilities.make_dir(bulk_out_dir)

        # symlink the enc file
        enc_name     = args.enc_ukb.split('/')[-1]
        ukb_enc_dest = os.path.join(bulk_out_dir, enc_name)
        utilities.make_symlink(args.enc_ukb, ukb_enc_dest)


        # build command for slurm submission
        convert_cmd  = 'cd {}'.format(bulk_out_dir)
        download_cmd = [args.ukb_conv, './{}'.format(enc_name), 'bulk', '-s{}'.format(bulk_id), '-o./{}'.format(bulk_id)]
        submit_cmd   = '{}\n{}'.format(convert_cmd, ' '.join(download_cmd))
        print(submit_cmd)

        write_file = os.path.join(args.slurm_dir, 'make_bulk_sublist_{}.bash'.format(bulk_id))
        utilities.write_cmd(cmd=submit_cmd, write_file=write_file)


if __name__ == "__main__":
    main()




