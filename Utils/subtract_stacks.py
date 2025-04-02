"""
subtract_stacks.py

Author: Thomas Plunkett

Date: 03/12/2024

Purpose: Subtract stacked images from each by first didiving into 4 tiles and then using OIS package
with Adaptive Bramich algorithm.
"""
# Import necessary packages 
import glob
import numpy as np
import pandas as pd
import argparse
import os
import astroalign as aa
from astropy.io import fits
from astropy.time import Time
from astropy import wcs
from astropy.wcs.utils import pixel_to_skycoord
import sep
from ois import optimal_system, eval_adpative_kernel
from scipy import ndimage, stats
import time
import datetime
import warnings
warnings.filterwarnings("ignore")

def register_tiles(im_list, ref_im, verbose):
    """
    A function to register images to the same pixel coordinate grid using AstroAlign.
    """
    if verbose:
        print('\n ------------------------------------------------------------------------ \n')
        print('\n COMMENCING REGISTRATION AND TILING! \n')
        print('\n ------------------------------------------------------------------------ \n')
        
    # Load in the reference frame data
    ref_data = fits.getdata(ref_im).astype(float)
    ref_header = fits.getheader(ref_im)
    ref_wcs = wcs.WCS(ref_header)
    ref_mask = ref_data > 30000
    ref_mask = ndimage.binary_dilation(ref_mask, iterations=10)
    masked_ref = np.ma.masked_array(ref_data, ref_mask)
    
    # Iterate through frames in list and register
    for im in im_list:
        if verbose:
            print('Registering frame: {:} \n'.format(im))
        
        # Get the header information and data from current image
        im_header = fits.getheader(im)
        source = fits.getdata(im).astype(float)
        source_wcs = wcs.WCS(im_header)
        source_mask = source > 30000
        source_mask = ndimage.binary_dilation(source_mask, iterations=10)
        masked_source = np.ma.masked_array(source, source_mask)
        
        # Check if WCS is present and remove, as now invalid.
        if source_wcs.has_celestial:
            try:
                del im_header[7:20]
            except:
                pass
        
        # If we aren't working on the ref im, use astroalign. Otherwise, just crop to subframe.
        if im != ref_im:
            try:
                count = 1
                registered_image, footprint = aa.register(masked_source, masked_ref,\
                                                          detection_sigma = 10, max_control_points = 200)
                
                # Divide image into 4 tiles to make computation faster and more robust
                M = registered_image.shape[0]//2
                N = registered_image.shape[1]//2
                for x in range(0,registered_image.shape[0],M):
                    for y in range(0,registered_image.shape[1],N):
                        reg_im = registered_image[x:x+M,y:y+N]
                        new_wcs = ref_wcs[x:x+M,y:y+N]
                        im_header.update(new_wcs.to_header(relax=True))
                        new_im = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_T{:}_reg.fits'.format(count)))
                        fits.writeto(new_im, reg_im, im_header, overwrite = True)
                        count += 1
            except:
                print('Failed to register frame: {:}... Moving on! \n'.format(im))
        else:
            count = 1
            M = source.shape[0]//2
            N = source.shape[1]//2
            # Divide image into 4 tiles to make computation faster and more robust
            for x in range(0, source.shape[0], M):
                for y in range(0, source.shape[1], N):
                    source_new = source[x:x+M,y:y+N]
                    new_wcs = ref_wcs[x:x+M,y:y+N]
                    ref_header.update(new_wcs.to_header(relax=True))
                    new_im = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_T{:}_ref.fits'.format(count)))
                    fits.writeto(new_im, source_new, ref_header, overwrite = True)
                    count += 1
                                  
            
def perform_subtile(reg_list, ref_im, kern_size, verbose):
    """
    Function to perform the image subtraction using the Adaptive Bramich algorithm.
    Saves to a new .fits file with _diff.fits at the end. 
    
    params:
    
    reg_list - The list of paths to registered images
    ref_im - The registered ref_im
    
    """
    # Get the reference image data and mask saturated pixels
    ref_data = fits.getdata(ref_im)
    ref_data = np.array(ref_data, dtype = "float")
    
    # Create mask - find the mode and if that value occurs more than 50 times, probably the image border to trim.
    # Also, find any pixel greater than 30000 counts and mask around that using binary dilation. 
    ref_mode, ref_count = stats.mode(ref_data, axis = None)
    ref_sat = ref_data > 30000
    ref_sat = ndimage.binary_dilation(ref_sat, iterations=10)
    if ref_count > 50:
        ref_bound = (ref_data == ref_mode[0])
        ref_mask = ref_bound | ref_sat
    else:
        ref_mask = ref_sat
    
    # Calculate the background with SEP package, faster than native method in OIS package and closer to SExtractor method
    ref_bkg = sep.Background(ref_data, mask = ref_mask, bw = 100, bh = 100, fw = 10, fh = 10).back()
    masked_ref = np.ma.masked_array(ref_data - ref_bkg, ref_mask)
    
    # Iterate through the list
    for im in reg_list:
        if im != ref_im:
            if verbose:
                print('Working on frame: {:} \n'.format(im))
                
            # Get target image data and mask saturated pixels
            im_data = fits.getdata(im).astype('float')
            im_header = fits.getheader(im)
            im_mode, im_count = stats.mode(im_data, axis = None)
            im_sat = im_data > 30000
            im_sat = ndimage.binary_dilation(im_sat, iterations=10)
            if im_count > 10:
                im_bound = (im_data == im_mode[0])
                im_mask = im_bound | im_sat
            else:
                im_mask = im_sat
                
            background = sep.Background(im_data, mask = im_mask, bw = 100, bh = 100, fw = 10, fh = 10).back()
            masked_im = np.ma.masked_array(im_data - background, im_mask)
            
            # Subtraction time!
            diff_image, optimal_image, kernel, _ = optimal_system(masked_im, masked_ref,\
                                                               kernelshape = (kern_size,kern_size), bkgdegree = None,\
                                                               method = 'AdaptiveBramich', poly_degree = 2)
            
            # Write kernel size to header
            im_header['KERNEL'] = kern_size
            
            # Mask bad pixels and save!
            diff_image.data[diff_image.mask] = 0
            fits.writeto(str(im).replace('.fits', '_diff.fits'), diff_image.data, im_header)
            fits.writeto(str(im).replace('.fits', '_bkg.fits'), background)
                                  
            
if __name__ == '__main__':
    # Set-up the parser
    parser = argparse.ArgumentParser(description='Perform DIA analysis on images in folder')
    parser.add_argument('Path', help='The path to folder containing .fits files.')
    parser.add_argument('Ref', help='The file name of the reference image to use.')
    parser.add_argument('Size', help='The kernel size (use 7 for default).')
    parser.add_argument('Verbose', type=str, help='Want to know whats happening? (y/n)')
    pargs = parser.parse_args()
    
    # Safety
    target_dir = os.path.abspath(pargs.Path)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)
    
    # Check and make a new directory for DIA data
    dia_dir = os.path.join(target_dir, 'DIA')
    if not os.path.isdir(dia_dir):
        os.mkdir(dia_dir)
    
    # Find the stacks
    stacks = glob.glob(os.path.join(target_dir, '*_Stack_*.fits'))
    
    # Convert to 64 bit and trim the edges
    for im in stacks:
        if 'weight' not in im:
            data = fits.getdata(im).astype(float)
            new_data = data[100:(2048 - 100), 100:(2048 - 100)]
            head = fits.getheader(im)
            im_wcs = wcs.WCS(head)
            head['BITPIX'] = '-64'
            head['NAXIS1'] = new_data.shape[0]
            head['NAXIS2'] = new_data.shape[1]
            new_wcs = im_wcs[100:(2048 - 100), 100:(2048 - 100)]
            head.update(new_wcs.to_header(relax=True))
            new_file = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_trim.fits'))
            fits.writeto(new_file, new_data, head, overwrite = True)
        
    # Find the stacks
    stacks = glob.glob(os.path.join(dia_dir, '*_Stack_*.fits'))
    
    # Specify the reference and then register the tiles
    ref_name = os.path.basename(pargs.Ref.replace('.fits', '_trim.fits'))
    ref = os.path.join(dia_dir, ref_name)
    register_tiles(stacks, ref, True)
    
    # Now do the difference imaging
    if pargs.Verbose == 'y' or pargs.Verbose == 'Y':
        print('\n --------------------------------------------------------------------------------- \n')
        print('\n COMMENCING DIFFERENCE IMAGING ON TILES... GO GET A COFFEE, THIS MAY TAKE A WHILE! \n')
        print('\n --------------------------------------------------------------------------------- \n')
        
    for i in range(1, 5):
        ref_im = glob.glob(os.path.join(dia_dir, '*_T{:}_ref.fits'.format(i)))[0]
        reg_im = glob.glob(os.path.join(dia_dir, '*_T{:}_reg.fits'.format(i)))
        
        if pargs.Verbose == 'n' or pargs.Verbose == 'N':
            perform_subtile(reg_im, ref_im, int(pargs.Size), False)
        else:
            perform_subtile(reg_im, ref_im, int(pargs.Size), True)
            
    if pargs.Verbose == 'y' or pargs.Verbose == 'Y':
        print('\n --------------------------------------------------------------------------------- \n')
        print('\n ALL DONE! HAVE A GOOD DAY =) \n')
        print('\n --------------------------------------------------------------------------------- \n')
        
