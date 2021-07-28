#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#

import os
import datetime
import calendar
import argparse
import logging
import utilities.utilities as utilities

log = logging.getLogger('ukb')

def main():
    '''
    This script runs "ukbunpack" on raw UKB enc files. 
    The actual ukbunpack command is a simple one-liner, but the added functionality
    makes sure that files are symlinked/populated in the correct place in the project directory.

    See here for more info: https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide

    Parameters
    ----------
    enc_key_pair: string, repeatable, required
        Full path to the UKB raw encoding file and key file, separated by colon
        [PATH_TO_ENC]:[PATH_TO_KEY]
        e.g. /Users/me/ukb12345.enc:/Users/me/ukb12345.key

    out_dir: string, required
        Where to write symlink and unpack data

    slurm: flag, optional
        Optionally submit the script directly to a cluster using SLURM 

    Returns
    -------
    command: string
        Path to a text file with the unpacking command. You can run this however you'd like. 

    '''

    # logging/feedback
    logging.basicConfig(level=logging.DEBUG)
    date = datetime.date.today()
    utilities.configLogging()
    log.info('buckner-lab-ukb-pipeline')
    log.info('RUNNING: UKB-DECRYPT')

    parser = argparse.ArgumentParser()
    parser.add_argument('--enc_key_pair', dest='enc_key_pair', required=True, action='append', nargs='+', default=None)
    parser.add_argument('--out_dir', dest='out_dir', required=True)
    parser.add_argument('--slurm', '-s', dest='slurm', action='store_true', default=False)
    args = parser.parse_args()

    # paths where data will be symlinked/decrypted
    raw_dir   = os.path.join(args.out_dir, 'data/ukb/raw')

    # ukb utility to unpack encrypted data
    ukbunpack = os.path.join(args.out_dir, 'data/ukb/external/ukbunpack')

    # decrypt each input enc/key pairs
    # ------
    for ukb_enc in args.enc_key_pair:

        enc_file = ukb_enc[0].split(':')[0]
        key_file = ukb_enc[0].split(':')[1]

        # 
        enc_name = enc_file.split('/')[-1]
        enc_base = enc_name.replace('.enc', '')

        # create destination enc file
        dest_enc = os.path.join(args.out_dir, enc_name)

        # make sure user fed the write filetype as input
        if enc_name[-4:] != '.enc':
            sys.exit('\n\nERROR: The UKB *enc path in your config file does not end in *.enc: {}'.format(enc_file))
        else:
            utilities.make_symlink(src_path=enc_file, dest_path=dest_enc)

        # symlink destination key file
        key_name = key_file.split('/')[-1]
        dest_key = os.path.join(args.out_dir, key_name)
        utilities.make_symlink(src_path=key_file, dest_path=dest_key)

        # construct the actual decryption command
        decrypt_cmd = f"cd {args.out_dir}\n"
        decrypt_cmd = f"{decrypt_cmd}{ukbunpack} ./{enc_name} ./{key_name}"

        # write command to file
        write_file = os.path.join(args.out_dir, 'decrypt_{}'.format(enc_base))
        log.info(write_file)

        # submit to slurm cluster 
        if args.slurm == True:
            utilities.submit_slurm(
                utilities.write_slurm(slurm_file=write_file,
                                      partition=args.slurm_partition,
                                      cmd=decrypt_cmd,
                                      jobName=enc_base,
                                      stime='6:00:00',
                                      nthreads=2))
        else:
            utilities.write_cmd(cmd=decrypt_cmd, write_file='{}.bash'.format(write_file))

if __name__ == "__main__":
    main()
