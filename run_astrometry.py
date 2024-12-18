"""
run_astrometry.py

Author: Thomas Plunkett

Date: 20/10/23

Purpose:

Perform astrometry on all images within a folder and subfolders

"""
# Import necessary packages 
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import update, getdata, getheader
from astropy.time import Time
import numpy as np
import os
import argparse
import multiprocessing as multi
from multiprocessing import cpu_count, TimeoutError
import subprocess as sub
import time
from prose import FitsManager

def get_center(img):
    """
    A function to find the center of the image in RA, DEC coords
    
    params:
    img - The path to the image
    
    return:
    center - A list containing the RA and DEC of the center
    """

    header = getheader(img)
    ra, dec = header["OBJCTRA"], header["OBJCTDEC"]
    date_img = Time(header['DATE-OBS'], format='isot', scale='utc')
    center = SkyCoord(ra, dec, unit=[u.hourangle, u.deg], obstime=date_img).fk5
    center = [center.ra.value, center.dec.value]

    return center

def set_wcs(img, pix_size):
    """
    A function to set the World Coordinate System for an image.
    
    params:
    
    img - The path to the image 
    pix_size - The pixel scale in arcseconds
    
    """
    
    # This function is based on astrometry.net, let's start to print some information
    #nameproc = int(multi.current_process().name.split("-")[1]) - 1
    start = time.localtime(time.time())
    print('Processing', os.path.basename(img), 'at', time.asctime(start))

    pix_low, pix_high = str(pix_size * 0.99), str(pix_size * 1.01)
    
    center = get_center(img)
    print(center)
    command = '/usr/local/astrometry/bin/solve-field --fits -Op -l 300 -t 3 -5 2.0 -u app -z 2' + ' -L ' + pix_low + ' -H ' + pix_high + ' -3 ' + str(center[0]) + ' -4 ' + str(center[1]) + ' ' + img
    
    process = sub.Popen([command], shell=True)
    process.wait()

    # If we don't solve, rename to indicate a bad file
    if not os.path.exists(img.replace('fits', 'solved')):
        bad_img = img.replace('.fits', '.bad')
        os.rename(img, bad_img)
        to_rm = img.replace('.fits', '.axy')
        os.remove(to_rm)

    new_img = img.replace('.fits', '.new')
    os.rename(new_img, img)
    
    for suf in ('-indx.xyls', '.axy', '.corr', '.match', '.rdls', '.solved', '.wcs'):
        to_rm = img.replace('.fits', suf)
        os.remove(to_rm)

    return

def check4wcs(im):
    """
    A function to check if astrometry has already been done on image.
    
    params:
    im - The path to the image to check
    """
    source = fits.open(im)[0]
    wcs_s = WCS(source.header)
    
    if wcs_s.has_celestial:
        return True 
    else:
        return False
        

if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Perform astrometry on images in a directory')
    parser.add_argument('object_folder', type=str, help='The path to the object directory')
    parser.add_argument('telescope', type=str, help='The telescope that collected this data')
    parser.add_argument('redo', type=str, help='Redo astrometry, even if WCS is found? (y/n)')
    args = parser.parse_args()
    
    # Safety and checks
    target_dir = os.path.abspath(args.object_folder)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)
    
    # Start action
    fm = FitsManager(target_dir, depth=5)
    
    if args.telescope == 'H50' or 'Planewave 50cm' or 'Harlingten 50cm':
        pix_size = 0.8
    else:
        print('Unidentified telescope!')
        raise SystemExit(1)
    
    # Force redo of astrometry
    if args.redo == 'Y' or args.redo == 'y':
        for img in fm.all_images:
            try:
                set_wcs(img, pix_size)
            except:
                print('Unable to find solution for {:}... Moving on!'.format(img))
    
    # Otherwise, check for WCS and only perform astrometry on images without
    else:
        for img in fm.all_images:
            try:
                if check4wcs(img) == False:
                    set_wcs(img, pix_size)
                else:
                    print('Image {:} already has astrometry. Moving on!'.format(img))
            except:
                print('Unable to find solution for {:}... Moving on!'.format(img))
