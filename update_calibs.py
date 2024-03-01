"""
update_calibs.py

Author: Thomas Plunkett

Purpose:

Check the nights folder for new calibration frames and move to relevant folders. Then create master frames!
"""
# Import necessary packages
from prose import Image, FitsManager
from astropy.io import fits
import numpy as np
import argparse
import glob
import os
from datetime import datetime 
import shutil
from UTGO_Pipeline.Utils.make_master_calibs import *

# Set up the parser
parser = argparse.ArgumentParser(description='Retieve calibration frames, organise them and create master frames')
parser.add_argument('night_folder', type=str, help='The path to the folder containing the nights data')
parser.add_argument('verbose', type=str, help='Do you want to know what is happening? (y/n)')
args = parser.parse_args()

# Define the paths
target_dir = os.path.abspath(args.night_folder)
dir_df = pd.read_csv('directories.csv')
calib_dir = os.path.abspath(dir_df.cal_dir[0])
flat_dir = os.path.join(calib_dir, 'Flats')
dark_dir = os.path.join(calib_dir, 'Darks')
bias_dir = os.path.join(calib_dir, 'Bias')
end_dir = os.path.abspath(dir_df.end_dir[0])

# Checks and safety
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)
    
if not os.path.isdir(calib_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Start the action!
fm = FitsManager(target_path, depth = 0)

# ----------------------------------------------------------------------------------------------- 
# Part 1 - Bias frames
if args.verbose == 'y' or 'Y':
    print('Looking for bias frames...')
    
if len(fm.all_bias) >= 15:
    bias_date = get_date_str(fm.all_bias[0])
    bias_fldr = os.path.join(bias_dir, bias_date)
    
    if not os.path.isdir(bias_fldr):
        os.mkdir(bias_fldr)
    
    # Copy to bias folder
    for bias in fm.all_bias:
        shutil.copy2(bias, bias_fldr)
    
    if args.verbose == 'y' or 'Y':
        print('Bias frames found! Creating master bias using sigma clipping...')
    
    m_bias_name = 'Master_Bias_'+date_str+'_sc.fits'    
    m_bias_path = os.path.join(bias_fldr, m_bias_name)
    make_master_bias(fm.all_bias, 'sigmaclip', m_bias_path)
    
    # Remove old master and copy
    master_list = glob.glob(os.path.join(bias_dir, 'Master_Bias_*'))
        if len(master_list) != 0:
            for m in master_list:
                os.remove(m)
    shutil.copy2(m_bias_path, bias_dir)    
    
# -----------------------------------------------------------------------------------------------       
# Part 2 - Dark frames
if args.verbose == 'y' or 'Y':
    print('Looking for dark frames...')
    
df_dark = fm.files(type='dark', path = True)
dark_date = get_date_str(fm.darks[0])
exp_list = df_dark['exposure'].unique()

for exp in exp_list:
    dark_list = df_dark[df_dark.exposure == exp]['path']
    if len(dark_list) >= 10:
        exp_str = str(int(exp))+'s'
        exp_dir = os.path.join(dark_dir, exp_str)
        dark_fldr = os.path.join(exp_dir, dark_date)
        
        if not os.path.isdir(exp_dir):
            os.mkdir(exp_dir)
        if not os.path.isdir(dark_fldr):
            os.mkdir(dark_fldr)
      
        for dark in dark_list:
            # Copy to darks folder
            shutil.copy2(dark, dark_fldr)
        
        if args.verbose == 'y' or 'Y':
            print('Dark frames found! Creating '+exp_str+' master dark using sigma clipping...')
        
        m_dark_name = 'Master_Dark_'+date_str+'_'+exp_str+'_sc.fits'    
        m_dark_path = os.path.join(dark_fldr, m_dark_name)
        make_master_dark(fm.all_darks, 'sigmaclip', m_dark_path)
        
        # Remove old master and copy
        master_list = glob.glob(os.path.join(exp_dir, 'Master_Dark_*'))
        if len(master_list) != 0:
            for m in master_list:
                os.remove(m)
        
        shutil.copy2(m_dark_path, exp_dir)    

# -----------------------------------------------------------------------------------------------       
# Part 3 - Flat frames
if args.verbose == 'y' or 'Y':
    print('Looking for flat frames...')
    
df_flat = fm.files(type='flat', path = True)
fltr_list = df_flat['filter'].unique()

# Get date and look for old folder
flat_date = get_date_str(fm.all_flats[0])
flat_monyr = flat_date[2:8]
old_dir = os.path.join(flat_dir, str('01'+flat_monyr))

if not os.path.isdir(old_dir):
    os.mkdir(old_dir)
    
# Check for repeats
fm_old = FitsManager(old_dir, depth=0)
old_df = fm_old.files(type='flat', path=True)
old_fltr_list = old_df['filter'].unique()
    
# If any repeats, make a new folder by date of flats
if bool(set(fltr_list).intersection(old_fltr_list)):
    for fltr in fltr_list:
        df_fltr = df_flat[df_flat['filter'] == str(fltr)]
        exp_list = df_fltr['exposure'].unique()
    
        for exp in exp_list:
            flat_list = df_fltr[df_fltr.exposure == exp]['path']
            flat_fldr = os.path.join(flat_dir, flat_date)
        
            # Copy to flats folder
            if len(flat_list) >= 5:
                if args.verbose == 'y' or 'Y':
                    print('Flat frames in {:} found! Creating master flats using median...'.format(fltr))
                    
                for flat in flat_list:
                    shutil.copy2(flat, flat_fldr)
                
                m_flat_name = 'Master_Flat_'+date_str+'_'+str(fltr)+'_med.fits'    
                m_flat_path = os.path.join(flat_fldr, m_flat_name)
                make_master_dark(flat_list, 'median', m_flat_path)
    
# Otherwise, default to old directory and fill with remaining flats
else:
    for fltr in fltr_list:
        df_fltr = df_flat[df_flat['filter'] == str(fltr)]
        exp_list = df_fltr['exposure'].unique()
    
        for exp in exp_list:
            flat_list = df_fltr[df_fltr.exposure == exp]['path']
            flat_fldr = old_dir
        
            # Copy to flats folder
            if len(flat_list) >= 5:
                if args.verbose == 'y' or 'Y':
                    print('Flat frames in {:} found! Creating master flats using median...'.format(fltr))
                    
                for flat in flat_list:
                    shutil.copy2(flat, flat_fldr)
                
                m_flat_name = 'Master_Flat_'+date_str+'_'+str(fltr)+'_med.fits'    
                m_flat_path = os.path.join(flat_fldr, m_flat_name)
                make_master_flat(flat_list, 'median', m_flat_path)

# ----------------------------------------------------------------------------------------------- 
# Fin!

if args.verbose == 'y' or 'Y':
    print('Calibration frames have been updated!')    