"""
uploader.py

Author: Thomas Plunkett

Org: UTAS

Purpose: Uploads data to Arcsecond.io 

"""
# Import necessary packages 
import pandas as pd
import numpy as np
import os
from pathlib import Path
import argparse
import glob
from arcsecond import ArcsecondConfig, UploadContext, FileUploader

# Set up the parser
parser = argparse.ArgumentParser(description='Code to upload data to Arcsecond.io portal')
parser.add_argument('target_dir', type=str, help='The path to the directory containing reduced images')
parser.add_argument('dataset', type=str, help = 'The dataset to associate the target to (i.e Transit, uLens, etc)')
parser.add_argument('clean', type=str, help = 'Upload only the best files? (y/n)')
parser.add_argument('verbose', type=str, help='Do you want to know what is happening? (y/n)')
args = parser.parse_args()
    
# Safety and checks
target_dir = Path(args.target_dir)
summary_dir = os.path.join(target_dir, 'Nightly_Summaries')
dir_df = pd.read_csv('/home/obs/UTGO_Pipeline/directories.csv')
uuid = dir_df.tele_uuid[0]

if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Set up the uploader
if args.verbose == 'y' or args.verbose == 'Y':
    print('Setting up the uploader... One moment please. \n')
    
config = ArcsecondConfig() # it will read your config file.
context = UploadContext(config, input_dataset_uuid_or_name=str(args.dataset), input_telescope_uuid=str(uuid),\
                        is_raw=False, org_subdomain="ariel-survey")
context.validate() # important step to perform before uploading.

# Find the good files to upload
if args.verbose == 'y' or args.verbose == 'Y':
    print('Searching for files to upload... Please wait! \n')

if os.path.isdir(summary_dir) and args.clean == 'y' or args.clean == 'Y':
    if args.verbose == 'y' or args.verbose == 'Y':
        print('Cleaning up bad files... \n')
    
    # Find all nightly summaries for this target
    sum_list = glob.glob(os.path.join(summary_dir, 'Summary_*.csv'))
    good_list = []
    
    # Iterate through list, add good files to master list to upload
    for s in sum_list:
        print(s)
        good_df = pd.read_csv(s)
        
        # Create a sub of good images, discarding files with FWHM or N_stars outside of 3 sigma from median
        stars_flag = np.round(good_df.N_stars.median() - 3*good_df.N_stars.std(), 0)
        fwhm_flag = np.round(good_df['FWHM [pix]'].median() + 3*good_df['FWHM [pix]'].std(), 2)
        bkg_flag = 50 #adu/s
        sub = np.logical_and(good_df.N_stars >= stars_flag, good_df['FWHM [pix]'] <= fwhm_flag)
        sub = np.logical_and(sub, good_df['Bkg [adu/s]'] <= bkg_flag)
        good_sub = good_df[sub]
        
        if args.verbose == 'y' or args.verbose == 'Y':
            bad_list = good_df[~sub]['File'].to_list()
            print('{:} subpar files found! Not uploading: {:}'.format(len(bad_list), bad_list))
            
        good_list += good_sub['File'].to_list()
        print(good_list)

# If we don't care about quality, still only upload files with astrometry
else:
    proc_csv = os.path.join(target_dir, 'proc_files.csv')
    good_df = pd.read_csv(proc_csv)
    good_list = good_df['Files'].to_list()

# Begin the uploading
if args.verbose == 'y' or args.verbose == 'Y':
    print('Beginning upload! \n')
    
for file in good_list:
    try:
        uploader = FileUploader(context, target_dir, Path(file), is_raw = False, display_progress=False)
        status, substatus, error = uploader.upload_file()
    except:
        print('Unable to upload {:}... Moving on!'.format(file))
        
# Delete the log files
log_files = glob.glob('/home/obs/UTGO_Pipeline/arcsecond_*.log')
for log in log_files:
    os.remove(log)
    
if args.verbose == 'y' or args.verbose == 'Y':
    print('Upload completed! \n')
