"""
flipper.py

Author: Thomas Plunkett

Purpose: Script to seperate meridian flipped images and rotate images by 180 degrees if necessary.
If astrometry already performed and flipping completed, will force a re-run of astrometry. 
"""

# Import necessary packages 
from astropy.io import fits
from astropy.wcs import WCS
import astroalign as aa
import numpy as np
import argparse
import os
import shutil
import glob
import warnings
import subprocess as sub

def find_flip(source_pth, target_pth):
    """
    Function to find if an image is flipped compared to a reference.
    
    params:
    source_pth - The path to the image to check
    target_pth - The path to the image to compare to
    
    return:
    is_flipped - Is the image flipped compared to target? (True or False)
    """
    source = fits.open(source_pth)[0]
    target = fits.open(target_pth)[0]
    wcs_s = WCS(source.header).wcs
    wcs_t = WCS(target.header).wcs
    
    # If we have WCS keywords, use those!
    if wcs_s.ctype[0] != '' and wcs_t.ctype[0] != '':
        if wcs_s.cd[0,0]*wcs_t.cd[0,0] > 0:
            is_flipped = False
        elif wcs_s.cd[0,0]*wcs_t.cd[0,0] < 0: 
            is_flipped = True
        elif wcs_s.cd[0,0]*wcs_t.cd[0,0] == 0:
            is_flipped = None 
            
    # No WCS, then use astroalign      
    else:
        try:
            transf, (source_list, target_list) = aa.find_transform(source.data, target.data)
            is_flipped = transf.rotation > np.pi/2
        except:
            is_flipped = None
    
    return is_flipped

def separate_flipped(target_dir, verbose):
    """
    A function to separate images into two folders (left and right)
    depending upon if it has been meridian flipped.
    
    params:
    target_dir - The folder to inspect for fits files
    verbose - Flag to tell if information is wanted ('y/n')
    """
    # Get a list of ims from targets directory
    os.chdir(target_dir)
    im_list = glob.glob('*.fits')
    target = im_list[0]
    
    # Make subfolders to contain images if they don't exist
    if not os.path.isdir('Left'):
        os.mkdir('Left')
    if not os.path.isdir('Right'):
        os.mkdir('Right')
    
    # Iterate through list and compare if flipped or not to a reference image,
    # which is just taken as first image in list
    for im in im_list[1:]:
        if verbose == 'y':
            print('Working on {:}'.format(im))
        if find_flip(im, target) == True:
            dst = os.path.join('Right', os.path.basename(im))
            shutil.move(im, dst)
            
        elif find_flip(im, target) == False: 
            dst = os.path.join('Left', os.path.basename(im))
            shutil.move(im, dst)
            
        else:
            print('Unable to find transform for {:}... Moving on!'.format(im))
    
    dst = os.path.join('Left', os.path.basename(im))
    shutil.move(target, dst)
    
    print('All done! Images separated into two folders')

def perform_flip(im, outpath):
    """
    Function to perform flipping and then save to a new .fits file.
    Note: currently does not flip WCS, if one exists. You may need to re-run astrometry.
    
    params:
    im - The image to flip
    """
    # Open file and get WCS, if it exists
    hdu = fits.open(im)[0]
    im_data = hdu.data
    
    # Do the flips...
    hdu.data = np.flip(im_data)
        
    # Add comment to header
    hdu.header['Comment'] = 'Image rotated by 180 degrees by flipper.py'

    # Write the cutout to a new FITS file
    flip_filename = os.path.basename(im)
    flip_outpath = os.path.join(outpath, flip_filename)
    hdu.writeto(flip_outpath, overwrite=True)
    
if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Separate images into folders based on side of meridian and then flip one set to match')
    parser.add_argument('path2fits', type=str, help='The path to .fits files to check and flip.')
    parser.add_argument('astr', type=str, help='Do the images have astrometry? (y/n)')
    parser.add_argument('verbose', type=str, help='Do you want to see the action? (y/n)')
    args = parser.parse_args()

    # Safety and checks
    target_path = os.path.abspath(args.path2fits)
    all_path = os.path.join(target_path, 'All')
    if not os.path.isdir(target_path):
        print("This folder doesn't exist!")
        raise SystemExit(1)
        
    if args.verbose != 'y' and args.verbose != 'n':
        print("Please pass only y or n to verbose flag!")
        raise SystemExit(1)
    
    # Separate into distinct folders, 'left' and 'right'
    separate_flipped(target_path, args.verbose)
    
    # Change to directory
    os.chdir(target_path)
    
    # Copy files from one folder to another
    if os.path.isdir('Left'):
        shutil.copytree('Left', 'All')
    
    # Flip the images on the other side of meridian and save to the All folder
    if os.path.isdir('Right') and os.path.isdir('All'):
        r_list = glob.glob(os.path.join('Right', '*.fits'))
        for rim in r_list:
            if args.verbose == 'y':
                print('Flipping file {:}'.format(rim))
            perform_flip(rim, all_path)
    
    # If there was astrometry on the original images, we need to redo it
    if args.astr == 'y' or args.astr == 'Y':
        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + os.path.join(target_path, 'Right') + ' H50 y'
        process = sub.Popen([cmd], shell = True)
        process.wait()

        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + os.path.join(target_path, 'All') + ' H50 y'
        process = sub.Popen([cmd], shell = True)
        process.wait()

    
    
