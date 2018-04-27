from . import utils
from . import file_parsers
from pyteomics import mgf as MGF
from tqdm import tqdm
import argparse

import os
from pathlib import Path


def get_mgfd_from_csv(csv):
    mses = file_parsers.load_frag_csv(csv)
    ds = [mse.to_mgf_dict() for mse in mses]    
    return ds

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input',help='Input file or folder to convert to MGF. Must be .csv')
    args= parser.parse_args()
    isdir = False
    if os.path.isdir(args.input):
        isdir = True
        csvs = [f.path for f in os.scandir(args.input) if f.name.lower().endswith('.csv')]
    else:
        csvs = [args.input]
    
    mgfs = [get_mgfd_from_csv(csv) for csv in tqdm(csvs,desc="Converting Files")]
    if isdir:
        output_folder = f"{args.input}_MGF"
        try:
            os.mkdir(output_folder)
        except OSError:
            print(f"Warning {output_folder} already exists, will be overwritten...")
        os.chdir(output_folder)
    for csv,mgf in tqdm(list(zip(csvs,mgfs)),desc='Writing Files'): 
        fname_root = os.path.splitext(os.path.split(csv)[1])[0]
        with open(f"{fname_root}.mgf",'w') as fout:
            MGF.write(mgf,output=fout)

if __name__ == "__main__":
    main()

        


        
