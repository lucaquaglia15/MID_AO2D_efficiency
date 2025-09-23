from telnetlib import Telnet
import matplotlib.pyplot as plt
import array as arr
import numpy as np
import os
import sys
import argparse
import yaml
import ROOT
import subprocess
from os import path
from pathlib import Path


def download(inputCfg):
    outDir = Path(inputCfg["input"]["output_dir_name"])
    run_list = outDir / "run_list.txt"

    if run_list.is_file():
        os.system(f'rm "{run_list}"')
        print(run_list)
        
    fOut = open(outDir / "run_list.txt", "x")

    print("----- Download and save files in %s -----" % (outDir))
    for iRun in range(0, len(inputCfg["input"]["run_list"])):

        file_type = inputCfg["input"]["file_type"]
        run = inputCfg["input"]["run_list"][iRun]
        alien_path = inputCfg["input"]["alien_input_path"][iRun]
        
        #Create folder with run number
        print("Run number: ",run)
        os.system(f'mkdir -p "{outDir / str(run)}"')
                
        #Copy from alien
        local_file = outDir / str(run) / "AnalysisResults.root"
        subprocess.run(["alien_cp",f"alien://{alien_path}/AnalysisResults.root",f"file:{local_file}"])

        fOut.write("{}\n".format(run))
    fOut.close()

def merge(inputCfg):
    fInPath = inputCfg["input"]["file_path"]
    file_type = inputCfg["input"]["file_type"]
    os.system("mkdir -p {}/merged_files".format(fInPath))
    runs = inputCfg["input"]["run_list"]
    for run in runs:
        #print("mkdir -p {}/merged_files/{}".format(fInPath, run))
        os.system(f'mkdir -p {fInPath}/merged_files/{run}')
        #print("hadd {}/merged_files/{}/AnalysisResults.root {}/{}/*/AnalysisResults.root".format(fInPath, run, fInPath, run))
        os.system(f'hadd {fInPath}/merged_files/{run}/{file_type} {fInPath}/{run}/*/{file_type}')

### ### ###
def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--download", help="Download single files", action="store_true")
    parser.add_argument("--merge", help="Do the merging of the downloaded files", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.download:
        download(inputCfg)
    if args.merge:
        merge(inputCfg)

main()