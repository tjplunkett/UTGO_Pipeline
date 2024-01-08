"""
prep4pysis.py

Author: Thomas Plunkett & JB Beaulieu 

Purpose:

Prepare images for difference imaging through PySIS. Code fixes header 
keywords, flips images either side of the meridian and, if selected, re-performs
astrometry. This is necessary if original files had astrometry.

"""
# Import necessary packages 
import os
import glob
import sys
import argparse
from flipper import *
from astropy.io import fits
from astropy.io.fits import getheader, getdata, update, open, setval
from astropy import time
from astropy.coordinates import SkyCoord, EarthLocation

if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Prepare reduced images for PySIS processing')
    parser.add_argument('path2fits', type=str, help='The path to .fits files to check and flip.')
    parser.add_argument('prefix', type=str, help='Image prefix to rename to')
    parser.add_argument('astr', type=str, help='Has atrometry been performed on the images? (y/n)')
    parser.add_argument('verbose', type=str, help='Do you want to see the action? (y/n)')
    args = parser.parse_args()

    # Safety and checks
    target_path = os.path.abspath(args.path2fits)
    all_path = os.path.join(target_path, 'All')
    prefix = args.prefix
    
    if not os.path.isdir(target_path):
        print("This folder doesn't exist!")
        raise SystemExit(1)
        
    if args.verbose != 'y' and args.verbose != 'n':
        print("Please pass only y or n to verbose flag!")
        raise SystemExit(1)
        
    if len(args.prefix) != 9:
        print("Make sure the prefix has 9 characters (i.e GB220022R)!")
        raise SystemExit(1)
    
    # Get a list of images
    im_list = glob.glob(os.path.join(target_path, '*.fits'))
    
    # Correct mjd
    print('Commencing correction of MJD & exposure keywords...')
    for im in im_list:
        if args.verbose == 'y':
            print("working on ",im)
        try:
            hdr = fits.open(im,mode='update')
        except:
            print("Error loading FITS file: " + im)
            sys.excepthook(*sys.exc_info())
            continue
    
        if args.verbose == 'y':
            print("Before modifications:")
            print()
            print(repr(hdr[0].header))
            print()
        # playing with keywords
        if args.verbose == 'y':
            print("values read ",hdr[0].header['EXPTIME'],hdr[0].header['DATE-OBS'])
    
        exptime = time.TimeDelta(hdr[0].header['EXPTIME'], format='sec').jd   # EXPTIME read in sec, then given its JD value
        jd_obs = time.Time(hdr[0].header['DATE-OBS'], format='isot',scale='utc').jd  # The same transformation for MJD-OBS from the appropriate format
        jd_mid = time.Time(jd_obs + exptime / 2., format='jd')
        mjd_obs = time.Time(hdr[0].header['DATE-OBS'], format='isot',scale='utc').mjd
    
        if args.verbose == 'y':
            print('EXPTIME=', exptime, 'jd_mid=', jd_mid,"MJD=",mjd_obs)

        hdr[0].header['PEXP'] = hdr[0].header['EXPTIME']
        hdr[0].header['PJD_MID'] = (jd_mid.value, 'keyword updated by Belial')
        hdr[0].header['MJD-OBS'] = (mjd_obs)
    
        hdr.flush()
        
        if args.verbose == 'y':
            print('all done')
 
    # Begin renaming
    print('Commencing renaming of files to match conventions...')
    if len(im_list) != 0:
        count = 1
        for im in im_list:
            # Work out the correct string length
            if count < 10:
                new_name = prefix+'000'+str(count)+'.fits'
            elif count > 9 and count < 100:
                new_name = prefix+'00'+str(count)+'.fits'
            elif count > 99 and count < 1000:
                new_name = prefix+'0'+str(count)+'.fits'
            elif count > 999 and count < 10000:
                new_name = prefix+str(count)+'.fits'
            
            dst = os.path.join(target_path, new_name)
            
            if args.verbose == 'y':
                print('File: {:} -> {:}'.format(im, dst))
            
            # Rename and iterate 
            os.rename(im, dst)
            count = count + 1
            
            
    # Finally perform the flippin'
    print('Commencing separation and flipping of images before and after meridian...')
     
    im_list = glob.glob(os.path.join(target_path, '*.fits'))
    
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
        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + os.path.join(target_path, 'Right') + ' H50'
        process = sub.Popen([cmd], shell = True)
        process.wait()

        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + os.path.join(target_path, 'All') + ' H50'
        process = sub.Popen([cmd], shell = True)
        process.wait()  
        