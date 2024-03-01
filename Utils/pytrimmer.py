"""
pytrimmer.py

Author: Thomas Plunkett

Date: 17/07/23

Purpose: Trim .fits images to required bounds

Usage: python pytrimmer.py [path2fits] [x_min] [x_max] [y_min] [y_max] [path2output]

"""
# Import necessary packages 
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import argparse
import os
import warnings
from astropy.utils.exceptions import AstropyWarning

# Define the function to trim and save images
def save_trim(im, x_min, x_max, y_min, y_max, outpath):
    """
    Function to perform trimming and then save to a new .fits file
    
    params:
    im - The path to the .fits files to trim
    x_min - The minimum x bound (int)
    x_max - The maximum x bound (int)
    y_min - The minimum y bound (int)
    y_max - The maximum y bound (int)
    outpath - The output path for saved images
    """
    # Open file and get WCS, if it exists
    hdu = fits.open(im)[0]
    im_wcs = WCS(hdu.header, relax=True)
    im_data = hdu.data
    
    # Check the max bounds aren't outside the image
    if x_max > im_data.shape[1]+1 or y_max > im_data.shape[0]+1:
        print("Make sure the max. bounds are within the image!")
        raise SystemExit(1)
        
    # Trim image
    trim_data = im_data[(y_min-1): (y_max-1), (x_min-1): (x_max-1)]
    trim_wcs = im_wcs[(y_min-1): (y_max-1), (x_min-1): (x_max-1)]
    
    # Put the cutout image in the FITS HDU
    hdu.data = trim_data
    
    # Update the FITS header with the cutout WCS
    hdu.header.update(trim_wcs.to_header(relax=True))
    
    # Add comment to header
    hdu.header['Comment'] = 'Image trimmed from {:}, with (x_min, x_max) = ({:},{:}) and (y_min, y_max) = ({:},{:})'.format(im, x_min, x_max, y_min, y_max)

    # Write the cutout to a new FITS file
    trim_str = '_'+str(trim_data.shape[1])+'x'+str(trim_data.shape[0])+'_trim.fits'
    trim_filename = os.path.basename(im.replace('.fits', trim_str))
    trim_outpath = os.path.join(outpath, trim_filename)
    hdu.writeto(trim_outpath, overwrite=True)

if __name__ == '__main__':
    # Suppress the FixFitsWarning
    warnings.simplefilter('ignore', AstropyWarning)

    # Set up the parser
    parser = argparse.ArgumentParser(description='Trim images to required bounds')
    parser.add_argument('im_path', type=str, help='The path to .fits files.')
    parser.add_argument('x_min', type=int, help='The minimum x-pixel')
    parser.add_argument('x_max', type=int, help='The maximum x-pixel')
    parser.add_argument('y_min', type=int, help='The minimum y-pixel')
    parser.add_argument('y_max', type=int, help='The maximum y-pixel')
    parser.add_argument('out_path', type=str, help='The path to save the new file to')
    args = parser.parse_args()

    # Safety and checks
    target_path = os.path.abspath(args.im_path)
    output_path = os.path.abspath(args.out_path)
    file_ext = os.path.splitext(target_path)[1]
    
    if not os.path.isfile(target_path):
        print("This file doesn't exist!")
        raise SystemExit(1)
    if not os.path.isdir(output_path):
        print("This output folder doesn't exist!")
        raise SystemExit(1)
    if file_ext != '.fits':
        print("This is not a .fits file!")
        raise SystemExit(1)
    if args.x_min > args.x_max or args.y_min > args.y_max:
        print("Make sure the min. bounds are smaller than max. bounds")
        raise SystemExit(1)
    if args.x_min < 1 or args.y_min < 1:
        print("Make sure the min. bounds are within the image (i.e >= 1)")
        raise SystemExit(1)

    # Run the action
    save_trim(target_path, args.x_min, args.x_max, args.y_min, args.y_max, output_path)
    print('Trimming complete!')