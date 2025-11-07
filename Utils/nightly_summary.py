"""
nightly_summary.py

Author: Thomas Plunkett

Date: 14/03/2024

Purpose:

Defines the functions needed to produce nightly summaries of observations. This
includes a .csv of nightly image stats (FWHM, Bkg, N_stars), a stack of the best 5 images 
and photometry catalogue from this stacked image. 

"""
# Import necessary packages
import os
import sys
import pandas as pd
import numpy as np
from prose import Image, Sequence, blocks, FitsManager, Observation
import subprocess as sub
from photutils.background import Background2D, SExtractorBackground
from astropy.io import fits
from datetime import timedelta
from astropy.time import Time

# Hacky shit to make directly callable
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Photometry.calibrate_phot import *

def writer(lst, title):
    """
    A function to write a list of files to an ascii file to pass to SWarp.
    
    params:
    lst - The list of files to write to ascii
    title - The title of the output ascii file
    """
    file_name = title
    with open(file_name, 'w') as f:
        for path in lst:
            f.write(f"{path}\n")
    f.close()

def make_stack(target_dir, objct, best_ims, bkg, date):
    """
    Creates a stacked image from the best images from a night.
    
    params:
    target_dir - The target directory
    objct - The object ID (str)
    best_ims - A list of the best images
    date - The date of observation
    """
    print(best_ims)
    
    if len(best_ims) != 0:
        cwd = os.path.abspath('/home/obs/UTGO_Pipeline')
        swarp_file = os.path.join(cwd, os.path.join('Utils', 'default.swarp'))
        list_file = str(objct+'_Stack_'+date+'.ascii')
        
        # Check for exisiting summary directory
        if not os.path.isdir(os.path.join(target_dir, 'Nightly_Summaries')):
            os.mkdir(os.path.join(target_dir, 'Nightly_Summaries'))
        
        # Write an ascii file for SWARP to read
        list_file = os.path.join(os.path.join(target_dir, 'Nightly_Summaries'), list_file)
        writer(best_ims, list_file)

        output_file = list_file.replace('.ascii', '.fits')

        # Run SWarp
        command = '/usr/bin/SWarp @' + list_file + ' -c ' + swarp_file + ' -IMAGEOUT_NAME ' + output_file + ' -WEIGHTOUT_NAME ' + output_file.replace('.fits', '.weight.fits')
        process = sub.Popen(str(command), shell = True, executable="/bin/bash")
        process.wait()

        # Add constant?
        if os.path.isfile(output_file):
           im_data = fits.getdata(output_file)
           stack_header = fits.getheader(output_file)
           stack_data = im_data + np.full_like(im_data, bkg)
           fits.writeto(output_file, data = stack_data, header = stack_header, overwrite = True)

        return output_file
           
def make_summary(target_dir, file_list, date):
    """
    Make a statistics file with FWHM, Bkg, etc for each image from a night of observing.
    
    params:
    target_dir - The target directory
    file_list - The list of file paths of images
    date - The date of observations
    """
    # Initialise lists to contain info
    bkg = []
    jd = []
    airmass = []
    good_ims = []
    best_ims = []
    n_stars = []
    best_bkgs = []
    med_fwhm = 0.0
    med_bkg = 0.0

    # Safety
    if not os.path.isdir(os.path.join(target_dir, 'Nightly_Summaries')):
        os.mkdir(os.path.join(target_dir, 'Nightly_Summaries'))
    
    # Define method to retrieve stats
    fwhm_data = blocks.Get("fwhm")
    get_stars = blocks.Get("sources")
    
    NSTARS = Sequence([
        blocks.detection.PointSourceDetection(), # stars detection
        get_stars,
    ])
    
    PSF = Sequence([
        blocks.detection.PointSourceDetection(n=10),  # stars detection
        blocks.Cutouts(clean = True), # making stars cutouts
        blocks.MedianPSF(),            # building PSF
        blocks.psf.Moffat2D(),     # modeling PSF
        fwhm_data,
    ])

    # Get object
    obj = fits.getval(file_list[0], 'OBJECT')
    obj = obj.replace('/','_')
    obj = obj.replace(' ', '_')

    
    # Loop through the images, test if they are good and then grab stats
    if len(file_list) != 0:
        for i in range(0, len(file_list)):
            try:
                NSTARS.run(file_list[i], show_progress=False)
                PSF.run(file_list[i], show_progress=False)
                n_stars += [len(get_stars.sources[i])]
                good_ims += [file_list[i]]
            except:
                print('Discarding image {:}...'.format(file_list[i]))
                os.rename(file_list[i], file_list[i].replace('.fits', '.bad'))
    
        # Next we do the background estimation!
    
        # Set our sigma clipper to +/- 3 st. devs
        sigma_clip = SigmaClip(sigma=3.0)

        if len(good_ims) != 0:
            for im in good_ims:
                img = Image(im)
                data = img.data
                hdr = img.header
                exp = float(hdr['EXPTIME'])
                
                try:
                    jd += [float(hdr['JD'])]
                except:
                    jd += [Time(hdr['DATE-OBS'], format='isot', scale='utc').jd]
                    
                airmass += [float(hdr['AIRMASS'])]
                sex_bkg = Background2D(data, (64, 64), filter_size=(9, 9), sigma_clip=sigma_clip, bkg_estimator=SExtractorBackground(), exclude_percentile = 25.0)   
                bkg += [sex_bkg.background_median/exp]
        
            # If more than 5 images, take just 5. If less, use all.
            if len(good_ims) >= 5 and len(fwhm_data.fwhm) != 0:        
                for i in np.argsort(fwhm_data.fwhm)[:5]:
                    best_ims += [good_ims[i]]
                    best_bkgs += [bkg[i]*exp]
            else: 
                best_ims = good_ims
                best_bkgs = np.array(bkg)*exp
        
        # Safety, in case FWHM estimation fails.
            if len(fwhm_data.fwhm) != len(good_ims):
                fwhm_data.fwhm = [np.nan]*len(good_ims)

            med_fwhm = np.median(fwhm_data.fwhm)
            med_bkg = np.median(best_bkgs)

        else:
            fwhm_data.fwhm = []
            n_stars = []
        
        # Create dataframe and write to .CSV
        summary_df = pd.DataFrame({'File': good_ims,'JD': jd, 'N_stars': n_stars, 'FWHM [pix]': fwhm_data.fwhm, 'Bkg [adu/s]': bkg, 'Airmass': airmass})
        summary_df = summary_df.round(3)
        summary_df.to_csv(os.path.join(os.path.join(target_dir, 'Nightly_Summaries'), 'Summary_H50_{:}_{:}.csv'.format(obj, date)))
                      
    return best_ims, med_fwhm, med_bkg

def run_sex(stack_image, ap_size, output_path):
    """
    A function to run SExtractor on stacked image and calibrate to GAIA DR2 synthetic phot.

    params:
    stack_image - Path to the stacked image
    ap_size - The size of the aperture to use (maybe make this the median FWHM)?
    output_path - Put into the nightly summary folder
    """
    cwd = '/home/obs/UTGO_Pipeline/'
    config_path = os.path.join(os.path.join(cwd, 'Photometry'), 'config')
    sex_path = os.path.join(config_path, 'default.sex')
    param_path = os.path.join(config_path, 'default_Taz50.param')
    cat_name = os.path.basename(stack_image).replace('.fits', '.cat')
    cat_path = os.path.join(output_path, cat_name)
    fltr = fits.getval(stack_image, 'FILTER')
    command = '/home/obs/astromatic/bin/sex ' + stack_image + ' -c '+ sex_path + ' -PARAMETERS_NAME ' + param_path + ' -PHOT_APERTURES ' + str(ap_size) + ' -CATALOG_NAME ' + cat_path
    process = sub.Popen([command], shell=True)
    process.wait()
                
    # Calibrate to GAIA
    df = read_sex(cat_path, stack_image)
    final_df = calibrate_phot(stack_image, df, fltr, output_path, True)
    final_df = final_df.round(4)
                              
    #Save file 
    file_name = str(os.path.basename(stack_image).replace('.fits','')+'_phot.csv')
    final_df.to_csv(os.path.join(output_path, file_name))

if __name__ == '__main__':
    
    # Set up the parser
    parser = argparse.ArgumentParser(description='Create nightly summary files')
    parser.add_argument('path', help='The path to folder containing .fits files.')
    parser.add_argument('date', type=str, help='The date of observations to produce files for. Type (i.e): 20250219 or all')
    parser.add_argument('stacks',type=str, help='Do you want to recreate the stacks and photometry? (y/n)')
    args = parser.parse_args()

    # Safety and checks
    target_dir = os.path.abspath(args.path)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)

    outpath = os.path.join(target_dir, 'Nightly_Summaries')

    # Begin the action!
    fm = FitsManager(target_dir, depth = 0)
    obs_df = fm.files(path = True)
    
    # To get the dates that are consistent with our naming convention (local times) need to some wizardry with time deltas
    date_list = obs_df.date.unique()
    
    if str(args.date) == 'all':
        # Iterate through all unique nights of observing, figure out the real date (local time) and run summary production
        for date in date_list:
            real_date = (pd.to_datetime(date) + timedelta(1)).strftime('%d%m%Y')
            print('Working on night: {:}'.format(real_date))
            obs_sub = obs_df[obs_df.date == date]
            file_sub = obs_sub.path.to_numpy()
            best_imgs, ap, bkg = make_summary(target_dir, file_sub, real_date)
            o = fits.getval(best_imgs[0],'OBJECT')

            if args.stacks == 'y' or args.stacks == 'Y':
                if len(best_imgs) > 1 and ap != 0.0 and bkg != 0:
                    try:
                        stack_im = make_stack(target_dir, o, best_imgs, bkg, real_date)
                        run_sex(stack_im, ap, os.path.join(target_dir, 'Nightly_Summaries'))
                        fits.setval(stack_im, keyword='APER_RAD', value=float(ap))
                    except:
                        raise
                        print('Unable to make stack for night {:} \n'.format(real_date))
                else:
                    print('Not enough images to stack on night {:} \n'.format(real_date))

    else:
        # Match the input date with the 'fake' dates from Prose, then create a subset of files to pass to the nightly summary function
        real_date = pd.to_datetime(args.date).strftime('%d%m%Y')
        fake_date = pd.to_datetime(str(args.date)) - timedelta(1)
        fake_date_str = fake_date.strftime('%Y-%m-%d')
        obs_sub = obs_df[obs_df.date == fake_date_str]
        file_sub = obs_sub.path.to_numpy()
        best_imgs, ap, bkg = make_summary(target_dir, file_sub, real_date)
        o = fits.getval(best_imgs[0],'OBJECT')

        if args.stacks == 'y' or args.stacks == 'Y':
                if len(best_imgs) > 1 and ap != 0.0 and bkg != 0:
                    try:
                        stack_im = make_stack(target_dir, o, best_imgs, bkg, real_date)
                        run_sex(stack_im, ap, os.path.join(target_dir, 'Nightly_Summaries'))
                    except:
                        print('Unable to make stack for night {:} \n'.format(real_date))
                else:
                    print('Not enough images to stack on night {:} \n'.format(real_date))

    print('Nightly summaries have been produced!')
            
