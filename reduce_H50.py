"""
reduce_H50.py

Author: Thomas Plunkett

Purpose:

Performs calibration (i.e dark-subtraction and flat-fielding) and then saves the reduced frame to a new file
"""

# Import necessary packages 
from prose import Image, Sequence, blocks, FitsManager, Observation
from astropy.io import fits
from Utils.pytrimmer import *
import numpy as np
import argparse
import glob
import os

# Constants and arrays
dim_arr = []
filter_arr = []
calib_list = []

# Set up the parser
parser = argparse.ArgumentParser(description='Reduce images')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
args = parser.parse_args()

target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Begin the action!
fm = FitsManager(target_dir, depth=args.depth)

# Change all the telescope headers to be compatible
for lit in fm.all_images:
    fits.setval(lit, 'TELESCOP', value='Planewave 50cm')
for drk in fm.all_darks:
    fits.setval(drk, 'TELESCOP', value='Planewave 50cm')
for bias in fm.all_bias:
    fits.setval(bias, 'TELESCOP', value='Planewave 50cm')
for flt in fm.all_flats:
    fits.setval(flt, 'TELESCOP', value='Planewave 50cm')

fm = FitsManager(target_dir, depth=args.depth)
print(fm)

calib_list += fm.all_darks.tolist()
calib_list += fm.all_bias.tolist()
calib_list += fm.all_flats.tolist()

# Check the dimensions of the science image and trim calibration frames as necessary
for im in fm.all_images:
        hdu = fits.open(im)[0]
        size_x = hdu.header['NAXIS1'] 
        size_y = hdu.header['NAXIS2']

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
     
    # Create a new subfolder to contain trimmed images, then trim 
    # and save the calibration frames
        if not os.path.isdir(path2dim):
            os.mkdir(path2dim)
            for cal in calib_list:
                save_trim(cal, subf_x, max_x, subf_y, max_y, path2dim)
        elif len(glob.glob(os.path.join(path2dim, '*.fit*'))) == 0:
            for cal in calib_list:
                save_trim(cal, subf_x, max_x, subf_y, max_y, path2dim)

# Loop over different filters and image sizes 
for pth in dim_arr:
    dm_str = os.path.basename(pth)
    wdth, hght = dm_str.split('x')
    fm_cal = FitsManager(pth, depth=0)
    for flt in filter_arr:
        flat_list = fm_cal.files(filter=flt,type='flat',path=True)['path'].to_numpy()
        im_list = fm.files(filter=flt, width=int(wdth), height=int(hght), type = 'light', path=True)['path'].to_numpy()
    
        # Define calibration sequence from raw files to reduced
        Calibrate = Sequence([
        blocks.Calibration(darks=fm_cal.all_darks, bias=fm_cal.all_bias, flats=flat_list),
        blocks.CleanBadPixels(darks = fm_cal.all_darks),
        blocks.SaveReduced(destination=pth, overwrite=True),
        ])
        
        if len(im_list) != 0:
            Calibrate.run(im_list)
        
        # Create a list of reduced images, then changes the image type to 'light' from 'reduced' 
        fits_list = glob.glob(os.path.join(pth,'*_reduced.fits'))
        for re in fits_list:
            fits.setval(re, 'IMAGETYP', value='Light')

print('Calibration complete!')
