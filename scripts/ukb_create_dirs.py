#!/usr/bin/python
#
# Author: Kevin Anderson, kevin.anderson@fas.harvard.edu
#

import os
import argparse


def make_dir(create_path):
    '''Create a path if it doesnt exist, but also print feedback'''
    if not os.path.exists(create_path):
        os.mkdir(create_path)
        print('Directory created: {}'.format(create_path))
    else:
        print('Directory exists: {}'.format(create_path) )


def main():
    '''
    Create project directories
    :param root_dir: directory where data will be written
    :return:
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--root_dir', dest='root_dir', required=True)
    args = parser.parse_args()

    make_dir(args.root_dir)
    dir_list = ['data', 
                'slurm', 
                'data/ukb',  
                'data/ukb/external', 
                'data/ukb/bulk', 
                'data/ukb/raw',
                'data/ukb/genetic', 
                'data/ukb/genetic/genotyped', 
                'data/ukb/genetic/imputed']
    for create_path in dir_list:
        make_dir(os.path.join(args.root_dir, create_path))

if __name__ == "__main__":
    main()

