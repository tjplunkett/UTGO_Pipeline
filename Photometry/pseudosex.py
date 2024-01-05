"""
pseudosex.py

Author: Thomas Plunkett

Purpose: 

Performs aperture photometry on a series of images using either SExtractor or Prose. 
This package uses the PhotUtils aperture photometry package by default. 

SourceExtractor will still be superior, but the PhotUtils does not require linux! 
Can be used on Windows.

"""
# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager, Observation
from calibrate_phot import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os
from astropy.wcs.utils import pixel_to_skycoord
from astropy.io import fits
import subprocess as sub

# Set up the parser 
parser = argparse.ArgumentParser(description='Perform source extraction on images within specified directories and output a catalog')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('ap_size', type=float, help='The size of the desired aperture for photometry? (use 5 if matching SExtractor)')
parser.add_argument('sex', type=str, help='Use SExtractor instead? (y/n)')
parser.add_argument('WCS_flag', type=str, help='Do these files have a WCS? (y/n)')
args = parser.parse_args()

# Check if the given directory exists
target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Get some paths
code_path = os.getcwd()
config_path = os.path.join(code_path, 'config')
sex_path = os.path.join(config_path, 'default.sex')
param_path = os.path.join(config_path, 'default_Taz50.param')

#Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
    
# Begin the action!
fm = FitsManager(target_dir, depth=args.depth)
print(fm)
image_list = fm.to_pandas('select path from files')['path'].to_numpy()

# image 0 will be our reference. Sometimes acts weird on stacks...
try:
    ref = Image(fm.all_images[0])
except:
    ref = Image(image_list[0])

detection = Sequence([
    blocks.Trim(),
    blocks.detection.PointSourceDetection(),
    blocks.Cutouts(),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])

# Run this and then show us the stars
detection.run(ref, show_progress=True)
ref.show()
plt.show()

data = blocks.Get("fluxes", "errors", "exposure", "airmass", 'sky', 'stars_coords', 'skycoord', 'wcs', 'filter')

# Define the photometry steps
photometry = Sequence([
    *detection[0:-1],        # apply the same detection process to all images
    blocks.detection.LimitStars(min=3),   # discard images not featuring enough stars
    blocks.BalletCentroid(),  
    blocks.PhotutilsAperturePhotometry(apertures = [args.ap_size], scale=False), # aperture photometry
    blocks.Peaks(),

    # Retrieving data from images in a conveniant way
    data,
])

# Loop through all images within the folder and perform source extraction, then try calibrate to GAIA DR2
for j in range(0, len(image_list)):
    try:
        if args.sex == 'y' or args.sex == 'Y':  
            try:
                # Try running source extractor
                cat_name = os.path.basename(image_list[j]).replace('.fits', '.cat')
                cat_path = os.path.join(output_path, cat_name)
                fltr = fits.getval(image_list[j], 'FILTER')
                command = 'sex ' + image_list[j] + ' -c '+ sex_path + ' -PARAMETERS_NAME ' + param_path\
                + ' -CATALOG_NAME ' + cat_path
                process = sub.Popen([command], shell=True)
                process.wait()
                
                # Calibrate to GAIA
                df = read_sex(cat_name)
                final_df = calibrate_phot(image_list[j], df, fltr, output_path)
                
            except:
                print('SExtractor has failed! Check your installation or use [sex] = n for default photomotry.')
        
        # Otherwise we want to use the default Prose photometry
        else:    
            photometry.run(image_list[j])
            sky = np.array(data.sky[j])
            fluxes, fluxes_er = np.array(data.fluxes[j]), np.array(data.errors[j])
            airmass = np.array(data.airmass[j])
            if airmass == None:
                airmass = 0
            mag, mag_er = -2.5*np.log10(fluxes), (2.5/np.log(10))*(fluxes_er/fluxes)
            pos = np.array(data.stars_coords[j])
            WCS = data.wcs[j]
            radec = np.array(data.skycoord[j])
            fltr = str(data.filter[j][0])
            length = len(fluxes[0])

            if args.WCS_flag == 'y':
                # Calculate the sky coordinates using WCS
                eq_coords = pixel_to_skycoord(pos[:,0], pos[:,1], WCS)

                d = {'RA': eq_coords.ra, 'DEC': eq_coords.dec, 'MAG_APER': mag[0], 'MAGERR_APER':mag_er[0], 'BACKGROUND': [sky]*length, 'X_IMAGE': pos[:,0], 'Y_IMAGE': pos[:,1], 'AIRMASS':[airmass]*length}

                # Calibrate to GAIA
                df = pd.DataFrame(d)
                final_df = calibrate_phot(image_list[j], df, fltr, output_path)

            else:
                d = {'X_IMAGE': pos[:,0], 'Y_IMAGE': pos[:,1], 'MAG_APER':  mag[0], 'MAGERR_APER':mag_er[0], 'BACKGROUND': [sky]*length, 'AIRMASS':[airmass]*length}
                final_df = pd.DataFrame(d)

        #Save file 
        file_name = str(os.path.basename(image_list[j]).replace('.fits','')+'_phot.csv')
        final_df.to_csv(os.path.join(output_path, file_name))
        
    except:
        raise
        print('Unable to perform source extraction on {:}... Moving on!'.format(image_list[j]))

print('Source extraction completed! Zeropoints have been written to images.')
