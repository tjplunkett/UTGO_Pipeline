"""
autoReduce.py

Author: Thomas Plunkett 
Organisation: UTAS, Physics
Date: 8/03/2024

Purpose: 

Perform automatic reduction on a folder containing a night of data. Called by reductionBot.py each day, but can be
run 'on the fly' in order to reduce data quickly during observations.

"""
# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager, Observation
from astropy.io import fits
from Utils.pytrimmer import *
from Utils.nightly_summary import *
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
parser.add_argument('on_the_fly', type=str, help="Run in 'on-the-fly' mode? (y/n)")
parser.add_argument('verbose', type=str, help='Do you want information on what is going on? (y/n)')
args = parser.parse_args()

# Define the paths
dir_df = pd.read_csv('/home/obs/UTGO_Pipeline/directories.csv')
target_dir = os.path.abspath(args.night_folder)
calib_dir = os.path.abspath(dir_df.cal_dir[0])
flat_dir = os.path.join(calib_dir, 'Flats')
dark_dir = os.path.join(calib_dir, 'Darks')
bias_dir = os.path.join(calib_dir, 'Bias')
end_dir = os.path.abspath(dir_df.end_dir[0])
multimess_dir = os.path.join(end_dir, 'MultiMess')
ulens_dir = os.path.join(end_dir, 'uLens')
transit_dir = os.path.join(end_dir, 'Transits')
other_dir = os.path.join(end_dir, 'Other')

# Define some arrays
transit_names = ['TOI','WASP', 'CoRoT','GJ', 'HAT', 'K2', 'KELT', 'KEPLER', 'NGTS', 'Qatar', 'TRAPPIST']
ulens_names = ['OGLE', 'OB', 'MOA', 'MB', 'KMT', 'UKIRT']
multimess_names = ['GRB', 'SN', 'GW']

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
    
# Define necessary functions
def event_type(ob):
    """
    A function to determine the type of event observed (i.e uLens, transit, etc).
    
    params: 
    ob - A string with the object ID
    
    return:
    fin_dir - The directory to send files to
    """
    fin_dir = other_dir 
    
    for name in transit_names:
        if name in ob:
            fin_dir = transit_dir  
    
    for name in ulens_names:
        if name in ob:
            fin_dir = ulens_dir
    
    for name in multimess_names:
        if name in ob:
            fin_dir = multimess_dir
    
    return fin_dir

# --------------------------------------------------------

### SECTION 1 - Find the right calibration frames to use 

# Determine which flats to use
if args.verbose == 'y':
    print('FINDING CLOSEST FLAT TO USE... \n')
    
# Find images
fm = FitsManager(target_dir, depth=0)

# Get the filters used on this night
light_df = fm.files(type='light')
fltr_list = light_df['filter'].unique()

# Retireve date and format correctly
hdu = fits.open(fm.all_images[0])[0]  
date_str = hdu.header['LOCALTIM'].split(' ')[0]
if len(date_str) == 7:
    date_str = '0'+date_str
im_date = datetime.strptime(date_str, '%d/%m/%Y').date()

# Venture through the flat folders to find closest by date and a complete set
flat_folders = glob.glob(os.path.join(flat_dir, '*'))
closest_flat = ''

for fldr in flat_folders:
    flat_str = os.path.basename(fldr)
    flat_date = datetime.strptime(flat_str, '%d%m%Y').date()
    
    # If taken before image, check if the folder contains all flats
    if flat_date <= im_date:
        flat_fm = FitsManager(fldr, depth = 0)
        flat_df = flat_fm.files(type='flat')
        flat_fltrs = flat_df['filter'].unique()
 
        if bool(set(fltr_list).issubset(flat_fltrs)):
            if closest_flat == '':
                closest_f_date = flat_date
                closest_flat = fldr
            elif im_date - flat_date < im_date - closest_f_date:
                closest_f_date = flat_date
                closest_flat = fldr 

# Onto finding the correct darks            
if args.verbose == 'y':
    print('FINDING THE CLOSEST DARKS... \n')

# Find all unique exposure times, find corresponding darks to match
exp_list = light_df['exposure'].unique()
dark_fldr_list = []
dark_exp_list = []
    
# Proceed to find closest darks (must be greater than or equal to sci im exp)
exp_folders = glob.glob(os.path.join(dark_dir, '*'))

# At the exposure level
for light_exp in exp_list:
    closest_exp = ''
    for exp_fldr in exp_folders:
        exp_str = os.path.basename(exp_fldr)
        dark_exp = float(exp_str.replace('s',''))

        if dark_exp == light_exp:
            closest_exp = dark_exp
        elif dark_exp > light_exp:
            if closest_exp == '':
                closest_exp = dark_exp
            elif dark_exp - light_exp < closest_exp - light_exp:
                closest_exp = dark_exp
    
    dark_folders = glob.glob(os.path.join(os.path.join(dark_dir, str(int(closest_exp))+'s'), '*'))
    closest_dark = ''
    
    # At the date level
    for dark_fldr in dark_folders:
        dark_str = os.path.basename(dark_fldr)
        dark_date = datetime.strptime(dark_str, '%d%m%Y').date()
        
        if closest_dark == '':
            closest_d_date = dark_date
            closest_dark = dark_fldr
        elif im_date - dark_date < im_date - closest_d_date:
                closest_d_date = dark_date
                closest_dark = dark_fldr
    
    dark_exp_list += [closest_exp]
    dark_fldr_list += [closest_dark]
    
# Find the correct bias to use - the easiest step
bias_folders = glob.glob(os.path.join(bias_dir, '*'))
closest_bias = ''

for bias_fldr in bias_folders:
    bias_str = os.path.basename(bias_fldr)
    bias_date = datetime.strptime(bias_str, '%d%m%Y').date()
        
    if closest_bias == '':
        closest_b_date = bias_date
        closest_bias = bias_fldr
    elif im_date - dark_date < im_date - closest_b_date:
        closest_b_date = bias_date
        closest_bias = bias_fldr

# ----------------------------------------------------
    
### SECTION 2 - Trimming calibration frames
if args.verbose == 'y':
        print('COMMENCING THE CALIBRATION PHASE... \n')

all_darks = []
master_darks = []

master_bias = glob.glob(os.path.join(closest_bias, 'Master_Bias_*.fits'))

for drk_fldr in dark_fldr_list:
    all_darks += glob.glob(os.path.join(drk_fldr, 'DARK_*.fits'))
    master_darks += glob.glob(os.path.join(drk_fldr, 'Master_Dark_*.fits'))
    
flat_frames = glob.glob(os.path.join(closest_flat, 'FLAT_*.fits'))
master_flats = glob.glob(os.path.join(closest_flat, 'Master_Flat_*.fits'))
    
calib_list = []
dim_arr = []
filter_arr = []
object_arr = []

# Begin the action!

# Change all the telescope headers to be compatible
for lit in fm.all_images:
    fits.setval(lit, 'TELESCOP', value='Planewave 50cm')
    
fm = FitsManager(target_dir, depth=0)

calib_list += master_flats
calib_list += master_bias
calib_list += master_darks
calib_list += all_darks
calib_list += flat_frames

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

    #filt = hdu.header['FILTER']
    max_x = subf_x + size_x
    max_y = subf_y + size_y
    dim_str = str(size_x)+'x'+str(size_y)
    path2dim = os.path.join(target_dir, dim_str)

    # Check for new dimensions
    if path2dim not in dim_arr:
        dim_arr += [path2dim]

    # Check for new filter
    #if filt not in filter_arr:
        #filter_arr += [filt]

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

#--------------------------------------------------------------------------------------
### SECTION 3 - Calibration/Reduction

# Create folders for each object to eventually move to
if not os.path.isdir(multimess_dir):
    os.mkdir(multimess_dir)

if not os.path.isdir(ulens_dir):
    os.mkdir(ulens_dir)
    
if not os.path.isdir(transit_dir):
    os.mkdir(transit_dir)
    
if not os.path.isdir(other_dir):
    os.mkdir(other_dir)

for o in object_arr:
    final_dir = event_type(o)
    obj_path = os.path.join(final_dir, o.replace(' ','_'))
    if not os.path.isdir(obj_path):
        os.mkdir(obj_path)

# Loop over different filters and image sizes 
if args.verbose == 'y':
    print('LOOPING THROUGH EACH SET OF FILTERS AND IMAGE SIZES... \n')

for pth in dim_arr:
    dm_str = os.path.basename(pth)
    wdth, hght = dm_str.split('x')
    master_bias = glob.glob(os.path.join(pth, 'Master_Bias_*.fits')) 
    
    # Tie the light image exposure time to the closest dark exposures and iterate
    for im_exp, dark_exp in zip(exp_list, dark_exp_list):
        for flt in fltr_list:
            flat_list = []
            master_flat = []
            flat_list = glob.glob(os.path.join(pth, str('FLAT_'+flt+'_*')))
            master_flat = glob.glob(os.path.join(pth, str('Master_Flat*_'+flt+'_*')))
            
            im_list = fm.files(filter=flt, exposure=im_exp, tolerance = 0, width=int(wdth), height=int(hght), type = 'light', path=True)['path'].to_numpy()
            
            # Check for images and FINALLY do the calibration!
            if len(im_list) != 0:
                dark_list = []
                master_dark = []
                dark_list = glob.glob(os.path.join(pth, str('DARK_*'+str(int(dark_exp))+'.000s*.fits')))
                master_dark = glob.glob(os.path.join(pth, str('Master_Dark*'+str(int(dark_exp))+'s*.fits')))             

                # Define calibration sequence from raw files to reduced
                if args.on_the_fly == 'y' and len(flat_list) != 0 and len(dark_list)!= 0:
                    Calibrate = Sequence([
                    blocks.Calibration(darks=master_dark, bias=master_bias, flats=flat_list),
                    blocks.SaveReduced(destination=pth, overwrite=True),
                    ])
                elif len(flat_list) != 0 and len(dark_list)!= 0:
                    Calibrate = Sequence([
                    blocks.Calibration(darks=master_dark, bias=master_bias, flats=flat_list),
                    blocks.CleanBadPixels(darks = dark_list),
                    blocks.SaveReduced(destination=pth, overwrite=True),
                    ])

                
                if len(im_list) != 0 and len(flat_list) != 0 and len(dark_list) != 0:
                    Calibrate.run(im_list)
                else:
                    print('UNABLE TO CALIBRATE THE FOLLOWING IMAGES: {:}. MOVING ON! \n'.format(im_list))

                # Create a list of reduced images, then changes the image type to 'light' from 'reduced' 
                # Also, sort into object folders as we go
                fits_list = glob.glob(os.path.join(pth,'*_reduced.fits'))
                comment_str = 'Calibrated using {:}, {:} and {:} by autoReduce.py'.format(master_dark, master_bias, master_flat)

                for re in fits_list:
                    fits.setval(re, 'IMAGETYP', value='Light')
                    fits.setval(re, 'SWCREATE', value = '')
                    hdu = fits.open(re)[0]
                    if 'HISTORY' not in hdu.header:
                        fits.setval(re, 'HISTORY', value=comment_str)
                    ob = hdu.header['OBJECT']
                    if ' ' in ob:
                        fits.setval(re, 'OBJECT', value=ob.replace(' ', '_'))
                        ob = ob.replace(' ', '_')
                    ob_flt = hdu.header['FILTER']
                    final_dir = event_type(ob)
                    ob_path = os.path.join(final_dir, ob)

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

#--------------------------------------------------------------------------------------
### SECTION 4 - Nightly Summaries and Astrometry (optional)

if args.on_the_fly == 'n' or args.on_the_fly == 'N': 
    if args.verbose == 'y':
        print('\n PERFORMING ASTROMETRY AND PRODUCING NIGHTLY SUMMARIES... \n')  
       
    for o in object_arr:
        final_dir = event_type(o)
        obj_path = os.path.join(final_dir, o)
        cmd = 'python /home/obs/UTGO_Pipeline/run_astrometry.py ' + str(obj_path) + ' H50 n'
        process = subprocess.Popen([cmd], shell=True)
        process.wait()
        
        # For each new object, gather files of each dimension and filter, create summaries
        for dim in glob.glob(os.path.join(obj_path, '*x*')):
            fm_object = FitsManager(dim, depth = 1)
            object_df = fm_object.files(type = 'light', path = True)
            obj_fltrs = object_df['filter'].unique()
            
            # Loop through filters, check they aren't narrowband
            for fltr in obj_fltrs:
                if fltr not in ['Ha', 'Hb', 'OIII', 'SII', 'Clear']:
                    file_lst = []
                    proc_files = []
                    obj_fltr_fldr = os.path.join(dim, fltr)
                    all_files = object_df[object_df['filter'] == str(fltr)]['path'].to_numpy()
                
                    # Check for existing list of processed files, if doesnt exist create a new one. 
                    # Otherwise, we just append to previous list.
                    if os.path.isfile(os.path.join(obj_fltr_fldr, 'proc_files.csv')):
                        proc_files = pd.read_csv(os.path.join(obj_fltr_fldr, 'proc_files.csv'))['Files'].to_list()
                        for fl in all_files:
                            if fl not in proc_files:
                                file_lst += [fl]
                               
                        if len(file_lst) != 0:   
                            best_imgs, ap, bkg  = make_summary(obj_fltr_fldr, file_lst, date_str.replace('/',''))
                            if len(best_imgs) > 1 and ap != 0.0 and bkg != 0:
                                try:
                                    stack_im = make_stack(obj_fltr_fldr, o, best_imgs, bkg, date_str.replace('/',''))
                                    run_sex(stack_im, ap, os.path.join(obj_fltr_fldr, 'Nightly_Summaries'))
                                    fits.setval(stack_im, keyword='APER_RAD', value=float(ap))
                                except:
                                    raise
                                    print('\n\n UNABLE TO PRODUCE SUMMARY FOR {:} IN {:} ON NIGHT: {:}... MOVING ON!'.format(o,fltr, date_str))
                            else:
                                print('\n\n NOT ENOUGH IMAGES! ABORTING SUMMARY FOR {:} IN {:} ON NIGHT: {:}... MOVING ON!'.format(o,fltr, date_str))
                                
                            proc_files += file_lst
                
                    else:
                        best_imgs, ap, bkg  = make_summary(obj_fltr_fldr, all_files, date_str.replace('/',''))
                        if len(best_imgs) > 1 and ap != 0.0 and bkg != 0:
                            try:
                                stack_im = make_stack(obj_fltr_fldr, o, best_imgs, bkg, date_str.replace('/',''))
                                run_sex(stack_im, ap, os.path.join(obj_fltr_fldr, 'Nightly_Summaries'))
                                fits.setval(stack_im, keyword='APER_RAD', value=float(ap))
                            except:
                                print('\n\n UNABLE TO PRODUCE SUMMARY FOR {:} IN {:} ON NIGHT: {:}... MOVING ON!'.format(o,fltr, date_str))
                        else:
                            print('\n\n NOT ENOUGH IMAGES! ABORTING SUMMARY FOR {:} IN {:} ON NIGHT: {:}... MOVING ON!'.format(o,fltr, date_str))
                            
                        proc_files = all_files
                
                    # Write out the proc_files.csv for book-keepin'
                    new_proc_df = pd.DataFrame({'Files': proc_files})
                    new_proc_df.to_csv(os.path.join(obj_fltr_fldr, 'proc_files.csv'))

                    #upload_cmd = 'python /home/obs/UTGO_Pipeline/uploader.py ' + str(obj_fltr_fldr) + ' ' + os.path.basename(final_dir) + ' y n'
                    upload_cmd = 'python /home/obs/UTGO_Pipeline/uploader.py ' + str(obj_fltr_fldr) + ' ' + str(o) + ' y n'
                    try:
                        upload_process = subprocess.Popen([upload_cmd], shell = True)
                        upload_process.wait()
                    except:
                        print('\n\n UNABLE TO UPLOAD TO ARCSECOND.IO PORTAL... MOVING ON!')
                        
                                         
if args.verbose == 'y':
    print('AUTO-REDUCTION IS COMPLETE! HISTORY HAS BEEN ADDED TO FITS HEADER.')    

        
    
