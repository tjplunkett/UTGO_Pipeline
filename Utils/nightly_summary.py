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
import pandas as pd
import numpy as np
import os
from prose import Image, Sequence, blocks, FitsManager, Observation
import subprocess as sub
from photutils.background import Background2D, SExtractorBackground
from astropy.stats import SigmaClip
from Photometry.calibrate_phot import *
from astropy.io import fits

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
        blocks.detection.PointSourceDetection(n=20),  # stars detection
        blocks.Cutouts(), # making stars cutouts
        blocks.MedianPSF(),            # building PSF
        blocks.psf.Moffat2D(),     # modeling PSF
        fwhm_data,
    ])

    # Get object
    obj = fits.getval(file_list[0], 'OBJECT')
    
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
        summary_df = pd.DataFrame({'File': good_ims, 'N_stars': n_stars, 'FWHM [pix]': fwhm_data.fwhm, 'Bkg [adu/s]': bkg})
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
    df = read_sex(cat_path)
    final_df = calibrate_phot(stack_image, df, fltr, output_path)
    final_df = final_df.round(4)
                              
    #Save file 
    file_name = str(os.path.basename(stack_image).replace('.fits','')+'_phot.csv')
    final_df.to_csv(os.path.join(output_path, file_name))
