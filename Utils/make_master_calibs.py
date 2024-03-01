"""
make_master_cals.py

Author: Thomas Plunkett

Purpose:

Creates master bias, darks and flat frames for use with the Harlingten 50cm telescope.

"""

# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager
from ccdproc import combine
from astropy.io import fits
import astropy.units as u
from astropy.stats import mad_std
import argparse
import os
import numpy as np

def get_date_str(image):
    """
    A convenience function to retrieve the date from an image
    """
    date_local = fits.getval(image, 'LOCALTIM')
    date_str = date_local.split(' ')[0]
    if len(date_str) == 7:
        date_str = '0'+date_str
    date_str = date_str.replace('/','')
           
    return date_str
            
def inv_median(a):
    """
    Function to caclulate the inverse of the median
    """
    return 1 / np.median(a)

def make_master_bias(im_list, method, output_filename):
    """
    A function to combine flat frames into a master flat.
    
    params:
    
    im_list - The list of flat images to combine
    """
    if method == 'sigmaclip':
        combined_bias = combine(im_list, method='average', \
            sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, \
                                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
    
    elif method == 'median':
        combined_bias = combine(im_list, method='median', unit = u.adu)
    
    combined_bias.meta['combined'] = True
    combined_bias.write(output_filename)

def make_master_darks(im_list, method, output_filename):
    """
    A function to combine flat frames into a master flat.
    
    params:
    
    im_list - The list of flat images to combine
    """
    if method == 'sigmaclip':
        combined_dark = combine(im_list, method='average',\
            sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, \
                                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
    
    elif method == 'median':
        combined_dark = combine(im_list, method='median', unit=u.adu)
    
    combined_dark.meta['combined'] = True
    combined_dark.write(output_filename)
    
    
def make_master_flat(im_list, method, output_filename):
    """
    A function to combine flat frames into a master flat.
    
    params:
    
    im_list - The list of flat images to combine
    """
    if method == 'sigmaclip':
        combined_flat = combine(im_list, method='average', scale=inv_median, sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=1.5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
        
    elif method == 'median':
        combined_flat = combine(im_list, method='median', scale=inv_median, unit=u.adu)
    
    combined_flat.meta['combined'] = True
    combined_flat.write(output_filename)
    

if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Create master calibration frames using ccdproc.')
    parser.add_argument('path', help='The path to folder containing .fits files.')
    parser.add_argument('method', help='The combination method (median or sigmaclip)')
    args = parser.parse_args()

    # Safety and checks
    target_dir = os.path.abspath(args.path)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)

    # Start action
    fm = FitsManager(target_dir)
    
    flat_df = fm.files(type='flat', path = True)
    flat_fltrs = flat_df['filter'].unique()
    
    if len(fm.all_bias) != 0:
        date_str = get_date_str(fm.all_bias[0])
        bias_name = 'Master_Bias_'+date_str+'.fits'
        bias_path = os.path.join(target_dir, bias_name)
        make_master_bias(fm.all_bias, args.method, bias_path)
        
    if len(fm.all_darks) != 0:
        date_str = get_date_str(fm.all_darks[0])
        dark_name = 'Master_Dark_'+date_str+'.fits'
        dark_path = os.path.join(target_dir, dark_name)
        make_master_darks(fm.all_darks, args.method, dark_path)
    
     
    if len(fm.all_flats) != 0:
        for fltr in flat_fltrs:
            print('Working on {:} flats...'.format(fltr))
            flat_list = flat_df[flat_df['filter'] == str(fltr)]['path']
            date_str = get_date_str(fm.all_flats[0]) 
            flat_name = 'Master_Flat_'+fltr+'_'+date_str+'.fits'
            flat_path = os.path.join(target_dir, flat_name)
            make_master_flat(flat_list, args.method, flat_path)

    print('Master frames have been created!')