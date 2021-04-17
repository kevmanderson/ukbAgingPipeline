
import os
import json
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', dest='config', required=True)
    opt = parser.parse_args()
    config_file = opt.config

    # read user configuration file
    # ------
    #config_file = '/gpfs/milgram/project/holmes/kma52/ukbAgingPipeline/config.json'
    with open(config_file, 'r') as f:
        config_json = json.load(f)[0]




        