"""
autoReduce.py

Author: Thomas Plunkett 
Organisation: UTAS, Physics
Date: 5/10/23

Purpose: 

Perform automatic reduction on a folder containing a night of data. Called by ReductionBot.py each day, but can be
run 'on the fly' in order to reduce data quickly during observations.

"""

# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager, Observation
from astropy.io import fits
from UTGO_Pipeline.Utils.pytrimmer import *
import numpy as np
import pandas as pd
import argparse
import glob
import os
from datetime import datetime 
import shutil
import subprocess

# Set up the parser
parser = argparse.ArgumentParser(description='Perform automatic reduction on images')
parser.add_argument('night_folder', type=str, help='The path to the folder containing the nights data')
parser.add_argument('astrometry', type=str, help='Perform astrometry during reduction? (y/n)')
parser.add_argument('verbose', type=str, help='Do you want information on what is going on? (y/n)')
args = parser.parse_args()

# Define the paths
dir_df = pd.read_csv('directories.csv')
target_dir = os.path.abspath(args.night_folder)
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
    
if not os.path.isdir(end_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# --------------------------------------------------------

### SECTION 1 - Find the right calibration frames to use 

# Read date from image header, determine which flats to use
if args.verbose == 'y':
    print('FINDING CLOSEST FLAT TO USE... \n')
          
fm = FitsManager(target_dir, depth=0)
hdu = fits.open(fm.all_images[0])[0]  
date_str = hdu.header['LOCALTIM'].split(' ')[0]

if len(date_str) == 7:
    date_str = '0'+date_str
im_date = datetime.strptime(date_str, '%d/%m/%Y').date()
flat_folders = glob.glob(os.path.join(flat_dir, '*'))
closest_flat = ''

for fldr in flat_folders:
    flat_str = os.path.basename(fldr)
    flat_date = datetime.strptime(flat_str, '%d%m%Y').date()
    if flat_date<im_date:
        if closest_flat == '':
            closest_date = flat_date
            closest_flat = fldr
        elif im_date - flat_date < im_date - closest_date:
            closest_date = flat_date
            closest_flat = fldr 

# Find max exposure time for night 
if args.verbose == 'y':
    print('FINDING THE CLOSEST DARK TO THE LONGEST EXPOSURE IMAGE... \n')
          
max_exp = 0
for im in fm.all_images:
    hdu = fits.open(im)[0]
    exp = float(hdu.header['EXPTIME'])
    if exp > max_exp:
        max_exp = exp 
    
# Proceed to find closest darks (must be greater than sci im exp)
dark_folders = glob.glob(os.path.join(dark_dir, '*'))
closest_dark = ''
closest_exp = 0

for drk_fldr in dark_folders:
    dark_str = os.path.basename(drk_fldr)
    dark_exp = float(dark_str.replace('s',''))
    
    if dark_exp>max_exp:
        if closest_dark == '':
            closest_dark = drk_fldr
            closest_exp = dark_exp
        elif dark_exp - max_exp < closest_exp - max_exp:
            closest_dark = drk_fldr
            closest_exp = dark_exp

if closest_dark == '':
    closest_dark = os.path.join(dark_dir, '600s')           

# ----------------------------------------------------
    
### SECTION 2 - Calibration
if args.verbose == 'y':
        print('COMMENCING THE CALIBRATION PHASE... \n')
          
# Commence le calibration
fm_bias = FitsManager(bias_dir, depth=0)
fm_dark = FitsManager(closest_dark, depth=0)
fm_flats = FitsManager(closest_flat, depth=0)
calib_list = []
dim_arr = []
filter_arr = []
object_arr = []

# Begin the action!

# Change all the telescope headers to be compatible
for lit in fm.all_images:
    fits.setval(lit, 'TELESCOP', value='Planewave 50cm')
for drk in fm_dark.all_darks:
    fits.setval(drk, 'TELESCOP', value='Planewave 50cm')
for bias in fm_bias.all_bias:
    fits.setval(bias, 'TELESCOP', value='Planewave 50cm')
for flt in fm_flats.all_flats:
    fits.setval(flt, 'TELESCOP', value='Planewave 50cm')

fm = FitsManager(target_dir, depth=0)

calib_list += fm_dark.all_darks.tolist()
calib_list += fm_bias.all_bias.tolist()
calib_list += fm_flats.all_flats.tolist()

# Check the dimensions of the science images and trim calibration frames as necessary
if args.verbose == 'y':
    print('CHECKING THE IMAGE SIZES AND TRIMMING CALIBRATION FRAMES... \n')
          
for im in fm.all_images:
        hdu = fits.open(im)[0]
        size_x = hdu.header['NAXIS1'] 
        size_y = hdu.header['NAXIS2']
        obj = hdu.header['OBJECT']

        if size_x == 2048:
            subf_x = 1
        else:
            subf_x = int(hdu.header['XORGSUBF'])+1
        if size_y == 2048:
            subf_y = 1
        else:
            subf_y = int(hdu.header['YORGSUBF'])+1 
        
        filt = hdu.header['FILTER']
        max_x = subf_x + size_x
        max_y = subf_y + size_y
        dim_str = str(size_x)+'x'+str(size_y)
        path2dim = os.path.join(target_dir, dim_str)
        
        # Check for new dimensions
        if path2dim not in dim_arr:
            dim_arr += [path2dim]
        
        # Check for new filter
        if filt not in filter_arr:
            filter_arr += [filt]
            
        # Check for new object
        if obj not in object_arr:
            object_arr += [obj]
     
    # Create a new subfolder to contain trimmed images, then trim 
    # and save the calibration frames
        if not os.path.isdir(path2dim):
            os.mkdir(path2dim)
            for cal in calib_list:
                save_trim(cal, subf_x, max_x, subf_y, max_y, path2dim)
        elif len(glob.glob(os.path.join(path2dim, '*.fit*'))) == 0:
            for cal in calib_list:
                save_trim(cal, subf_x, max_x, subf_y, max_y, path2dim)

# Create folders for each object to eventually move to
for o in object_arr:
    obj_path = os.path.join(end_dir, o)
    if not os.path.isdir(obj_path):
        os.mkdir(obj_path)

# Loop over different filters and image sizes 
if args.verbose == 'y':
    print('LOOPING THROUGH EACH SET OF FILTERS AND IMAGE SIZES... \n')

for pth in dim_arr:
    dm_str = os.path.basename(pth)
    wdth, hght = dm_str.split('x')                    
    fm_cal = FitsManager(pth, depth=0)
    
    for flt in filter_arr:
        flat_list = fm_cal.files(filter=flt,type='flat',path=True)['path'].to_numpy()
        im_list = fm.files(filter=flt, width=int(wdth), height=int(hght), type = 'light', path=True)['path'].to_numpy()
    
        # Define calibration sequence from raw files to reduced
        if len(flat_list) != 0:
            Calibrate = Sequence([
            blocks.Calibration(darks=fm_cal.all_darks, bias=fm_cal.all_bias, flats=flat_list),
            blocks.CleanBadPixels(darks = fm_cal.all_darks, flats = flat_list),
            blocks.SaveReduced(destination=pth, overwrite=True),
            ])
        
        if len(im_list) != 0 and len(flat_list) != 0:
            Calibrate.run(im_list)
        else:
            print('UNABLE TO CALIBRATE THE FOLLOWING IMAGES: {:}. CHECK IF FLATS ARE AVAILABLE! \n'.format(im_list))
        
        # Create a list of reduced images, then changes the image type to 'light' from 'reduced' 
        # Also, sort into object folders as we go
        fits_list = glob.glob(os.path.join(pth,'*_reduced.fits'))
        comment_str = 'Calibrated using {:} x {:} s darks, {:} bias frames and {:} {:} flats from {:} by autoReduce.py'.format(len(fm_cal.all_bias.tolist()), max_exp, len(fm_cal.all_bias.tolist()), len(flat_list), flt, closest_date)
        
        for re in fits_list:
            fits.setval(re, 'IMAGETYP', value='Light')
            hdu = fits.open(re)[0]
            if 'HISTORY' not in hdu.header:
                fits.setval(re, 'HISTORY', value=comment_str)
            ob = hdu.header['OBJECT']
            if ' ' in ob:
                fits.setval(re, 'OBJECT', value=ob.replace(' ', '_')
                ob = ob.replace(' ', '_')
            ob_flt = hdu.header['FILTER']
            ob_path = os.path.join(end_dir, ob)
            
            if ob == 'Entered Coordinates':
                ob_path = ob_path.replace(' ','_')
                date_path = os.path.join(ob_path, date_str.replace('/',''))
                obdim_path = os.path.join(date_path, dm_str)
                final_path = os.path.join(obdim_path, ob_flt)
                
                # Make directories if they dont exist
                if not os.path.isdir(date_path):
                    os.mkdir(date_path)
                if not os.path.isdir(obdim_path):
                    os.mkdir(obdim_path)
                if not os.path.isdir(final_path):
                    os.mkdir(final_path)  
            else:
                obdim_path = os.path.join(ob_path, dm_str)
                final_path = os.path.join(obdim_path, ob_flt)
                # Make directories if they dont exist
                if not os.path.isdir(obdim_path):
                    os.mkdir(obdim_path)
                if not os.path.isdir(final_path):
                    os.mkdir(final_path)
   
            shutil.copy(re, final_path)
    
    #delete the folder 
    if os.path.isdir(pth):
        shutil.rmtree(pth)

if args.astrometry == 'y' or args.astrometry == 'Y':        
    for o in object_arr:
        obj_path = os.path.join(end_dir, o)
        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + str(obj_path) + ' H50 n'
        process = subprocess.Popen([cmd], shell=True)
        process.wait()

if args.verbose == 'y':
    print('AUTO-REDUCTION IS COMPLETE! HISTORY HAS BEEN ADDED TO FITS HEADER.')    

        
    
