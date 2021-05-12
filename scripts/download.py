#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#
"""This module provides utilities for downloading and converting UKB data"""

import os
import sys
import wget

from utilities.utilities import make_symlink

def get_ukbutils(util_dir):
    '''Download UKB utilities from here: https://biobank.ndph.ox.ac.uk/ukb/download.cgi'''
    utils = ['ukbmd5', 'ukbconv', 'ukbunpack', 'ukbfetch', 'ukblink', 'gfetch']
    for ukb_util in utils:
        util_url = 'http://biobank.ndph.ox.ac.uk/ukb/util/{}'.format(ukb_util)
        dest_path= os.path.join(util_dir, ukb_util)
        if not os.path.exists(dest_url):
            print('Downloading util: {}'.format(os.path.join(util_dir, ukb_util)))
            wget.download(util_url, out=util_dir)
            os.chmod(dest_path, 755)

    # encoding file
    wget.download('biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb', out=util_dir)
    os.chmod(os.path.join(util_dir, 'encoding.ukb'), 755)
    

def decrypt_ukb_data(ukbunpack, write_dir, enc, key):

    # create destination enc file
    enc_name = enc.split('/')[-1]
    dest_enc = os.path.join(write_dir, enc_name)

    # make sure user fed the write filetype as input
    if enc_name[-4:] != '.enc':
        sys.exit('\n\nERROR: The UKB *enc path in your config file does not end in *.enc: {}'.format(enc))
    else:
        make_symlink(src_path=enc, dest_path=dest_enc)

    # create destination key file
    key_name = key.split('/')[-1]
    dest_key = os.path.join(write_dir, key_name)
    make_symlink(src_path=enc, dest_path=dest_key)

    decrypt_cmd = f"cd {write_dir}\n"
    decrypt_cmd = f"{decrypt_cmd}{ukbunpack} ./{enc_name} ./{key_name}"
    return decrypt_cmd


    # convert command
    print('Unpacking UKB Data')
    os.chdir(raw_dir)
    converted_data = './{}.enc_ukb'.format(enc_base)
    if not os.path.exists(converted_data) or opt.overwrite == True:
        if os.path.exists(converted_data):
            os.remove(converted_data)
        subprocess.call([os.path.join(config_json['repo_dir'], 'external/ukbunpack'), './{}.enc'.format(enc_base),
                         './{}'.format(key_name)])
        # subprocess.call(['/ukbtools/ukbunpack', './{}.enc'.format(enc_base), './{}'.format(key_name)])








