"""
fix_header.py 

Author: Thomas Plunkett

Purpose:

Bulk changing of fits keyword values or addition of new keywords/values to files
"""
# Import necessary packages 
from astropy.io import fits
import glob
import os
import argparse

# Main action! 
if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Bulk changing of fits keywords or addition of new keywords to files')
    parser.add_argument('path', type=str, help='The path to the folder containing .fits files')
    parser.add_argument('keyword', type=str, help='The keyword to target or create')
    parser.add_argument('value', help='The value to change or assign to the keyword')
    args = parser.parse_args()
    
    # Safety and checks
    target_dir = os.path.abspath(args.path)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)
    
    # Go to folder and make a list of fits files
    os.chdir(target_dir)
    fits_list = glob.glob('*.fit*')
    
    # If empty list, check for .FIT, which comes from 1.3 m
    if len(fits_list) == 0:
        fits_list = glob.glob('*.FIT')
    
    # Set the desired values in the header for each image
    for im in fits_list:
        fits.setval(im, keyword=args.keyword, value=args.value)

    print('All done!')