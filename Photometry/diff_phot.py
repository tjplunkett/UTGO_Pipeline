"""
diff_phot.py

Author: Thomas Plunkett
Organisation: UTAS
Date: 3/04/23

Purpose: This script performs differential photometry on calibrated images using the 'prose' python package
"""
# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager, Observation
import query_gaia
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os
import glob

# Set up the parser
parser = argparse.ArgumentParser(description='Perform differential photometry on images within specified directories and output LC and other plots')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('threshold', type=float, help='Source detection threshold (sigmas above background).')
parser.add_argument('obs_date', type=str, help='The date of observation (ddmmyyyy) or all')
args = parser.parse_args()

# Safety and checks
target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)

# Define path to nightly summaries
night_sum = os.path.join(target_dir, 'Nightly_Summaries')
proc_path = os.path.join(target_dir, 'proc_files.csv')

# Begin the action!
# ------------------------------------------------------------------------
# Part 1 - Find a good reference image, model PSF and discard bad images
fm = FitsManager(target_dir, depth=args.depth)
print(fm)

# Find a reference image from best seeing and discard bad images
print('Searching for reference image... Please wait!')

if os.path.isdir(night_sum):
    try:
        if args.obs_date == 'all' or args.obs_date == 'All':
            proc_df = pd.read_csv(proc_path)
            sum_file = glob.glob(os.path.join(night_sum, str('Summary_*.csv')))[0]
            sum_df = pd.read_csv(sum_file)
            good_ims = proc_df['Files']
            best_img = sum_df[sum_df['FWHM [pix]'] == sum_df['FWHM [pix]'].min()]['File'].to_list()[0]
            best_im = Image(best_img)
        else:
            sum_file = glob.glob(os.path.join(night_sum, str('Summary_*'+args.obs_date+'*')))[0]
            sum_df = pd.read_csv(sum_file)
            good_ims = sum_df['File'] 
            best_img = sum_df[sum_df['FWHM [pix]'] == sum_df['FWHM [pix]'].min()]['File'].to_list()[0]
            best_im = Image(best_img)
    except:
        raise
        print('No nightly summary found... Continuing manual search for best image...')

else:   
    fwhm_data = blocks.Get("fwhm")
    good_ims = []

    PSF = Sequence([
        blocks.detection.PointSourceDetection(n=50),  # stars detection
        blocks.Cutouts(), # making stars cutouts
        blocks.MedianPSF(),                 # building PSF
        blocks.psf.FastGaussian(),              # modeling PSF
        fwhm_data,
    ])

    for im in fm.all_images:
        try:
            print('Working on: {:}...'.format(im))
            PSF.run(im, show_progress=False)
            good_ims += [im]
        except:
            print('Discarding image {:}...'.format(im))

        best_img = good_ims[np.argmin(np.array(fwhm_data.fwhm))]
        best_im = Image(best_img)

print('Best image for night {:}: {:}'.format(args.obs_date, best_img))

detection = Sequence([
    blocks.detection.PointSourceDetection(threshold = args.threshold),  # stars detection
    blocks.Cutouts(clean = False), # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
    ])

detection.run(best_im, show_progress = False)    
best_im.show(stars=True)
plt.show()

# ----------------------------------------------------------------------------------------
# Part 2 - Perform the photometry on desired target

# Choose the target
target_no = int(input('Which is the target?'))

# Check if photometry already exists
if os.path.isfile(os.path.join(output_path, 'data.phot')):
    obs = Observation(os.path.join(output_path, 'data.phot'))
else:
   # Define the photometry steps
    photometry = Sequence([
        *detection[0:-1],  # apply the same calibration to all images
        blocks.psf.Moffat2D(reference=best_im),   # providing a reference improve the PSF optimisation
        blocks.detection.LimitStars(min=20),   # discard images not featuring enough stars
        blocks.Twirl(best_im.stars_coords, n=20),       # compute image transformation

       # set stars to the reference ones and apply the inverse
       # transformation previously found to match the ones in the image
        blocks.Set(stars_coords=best_im.stars_coords),
        blocks.AffineTransform(data=False, inverse=True),

        blocks.BalletCentroid(),              # stars centroiding
        blocks.PhotutilsAperturePhotometry(), # aperture photometry
        blocks.Peaks(),

        # Retrieving data from images in a conveniant way
        blocks.XArray(
            ("time", "jd_utc"),
            ("time", "bjd_tdb"),
            ("time", "flip"),
            ("time", "fwhm"),
            ("time", "fwhmx"),
            ("time", "fwhmy"),
            ("time", "dx"),
            ("time", "dy"),
            ("time", "airmass"),
            ("time", "exposure"),
            ("time", "path"),
            ("time", "sky"),
            (("time", "apertures", "star"), "fluxes"),
            (("time", "apertures", "star"), "errors"),
            (("time", "apertures", "star"), "apertures_area"),
            (("time", "apertures", "star"), "apertures_radii"),
            (("time", "apertures"), "apertures_area"),
            (("time", "apertures"), "apertures_radii"),
            ("time", "annulus_rin"),
            ("time", "annulus_rout"),
            ("time", "annulus_area"),
            (("time", "star"), "peaks"),
            name="xarray"
        ),

       # Stack image
        blocks.AffineTransform(stars=False, data=True),
        blocks.Stack(best_im, name="stack"),
    ])

    # Run the photometry on all images, then convert this to an obs object
    photometry.run(good_ims)
    obs = Observation(photometry.xarray.to_observation(photometry.stack.stack, sequence=photometry))

# Set star of interest
obs.target = target_no
    
# Choose the method to use to find comparison stars
if input('Would you like to query Gaia to find comparison star? (y/n)').lower() == 'y':
    try:
        print('Querying Gaia... This may take a moment')
        comp_idx = query_gaia.find_comp_bycol2(best_im, target_no)
        obs.diff(comp_idx)
        obs.plot_summary()
        plt.savefig(os.path.join(output_path, 'DR2_Summary'))
        obs.show_stars(size=5)
        plt.savefig(os.path.join(output_path, 'DR2_CompStars'))
        #obs.plot_comps_lcs()
        #plt.savefig(os.path.join(output_path, 'DR2_CompStarLC'))
        #obs.plot_precision()
        #plt.savefig(os.path.join(output_path, 'DR2_Precision'))
        obs.plot_psf_model()
        plt.savefig(os.path.join(output_path, 'PSF_Model'))
        plt.show()
    except:
        print('Problem with querying gaia... Using Broeg algorithm instead!')
        obs.broeg2005()
        obs.plot_summary()
        plt.savefig(os.path.join(output_path, 'Broeg_Summary'))
        obs.show_stars(size=5)
        plt.savefig(os.path.join(output_path, 'Broeg_CompStars'))
        #obs.plot_comps_lcs()
        #plt.savefig(os.path.join(output_path, 'Broeg_CompStarLC'))
        #obs.plot_precision()
        #plt.savefig(os.path.join(output_path, 'Broeg_Precision'))
        obs.plot_psf_model()
        plt.savefig(os.path.join(output_path, 'PSF_Model'))
        plt.show()

else: 
    print('Using Broeg algorithm instead... Please wait')
    obs.broeg2005()
    obs.plot_summary()
    plt.savefig(os.path.join(output_path, 'Broeg_Summary'))
    obs.show_stars(size=5)
    plt.savefig(os.path.join(output_path, 'Broeg_CompStars'))
    #obs.plot_comps_lcs()
    #plt.savefig(os.path.join(output_path, 'Broeg_CompStarLC'))
    #obs.plot_precision()
    #plt.savefig(os.path.join(output_path, 'Broeg_Precision'))
    obs.plot_psf_model()
    plt.savefig(os.path.join(output_path, 'PSF_Model'))
    plt.show()

# Save the results to a csv
obs.to_csv(os.path.join(output_path, 'diffLC.csv'), sep=',')
obs.save(os.path.join(output_path, 'data.phot'))
