"""
diffphot_pipe.py

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

# Set up the parser
parser = argparse.ArgumentParser(description='Perform differential photometry on images within specified directories and output LC and other plots')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
args = parser.parse_args()

# Safety and checks
target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Begin the action!
fm = FitsManager(target_dir, depth=args.depth)
print(fm)

# image 0 will be our reference
ref = Image(fm.all_images[0])

# Define 'calibration' steps, though images should already be calibrated
calibration = Sequence([
    blocks.SegmentedPeaks(min_area=1.5, minor_length=1),  # stars detection
    blocks.Cutouts(clean=False),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])

# Run this and then show us the stars to pick the target
calibration.run(ref, show_progress=True)
ref.show()
plt.show()

# Choose the target
target_no = int(input('Which is the target?'))

# Define the photometry steps
photometry = Sequence([
    *calibration[0:-1],                   # apply the same calibration to all images
    blocks.psf.Moffat2D(reference=ref),   # providing a reference improve the PSF optimisation
    blocks.detection.LimitStars(min=3),   # discard images not featuring enough stars
    blocks.Twirl(ref.stars_coords),       # compute image transformation

    # set stars to the reference ones and apply the inverse
    # transformation previously found to match the ones in the image
    blocks.Set(stars_coords=ref.stars_coords),
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
    blocks.Stack(ref, name="stack"),
])

# Run the photometry on all images, then convert this to an obs object
photometry.run(fm.all_images)
obs = Observation(photometry.xarray.to_observation(photometry.stack.stack, sequence=photometry))
obs.target = target_no

# Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
    
# Choose the method to use to find comparison stars
if input('Would you like to query Gaia to find comparison star? (y/n)').lower() == 'y':
    try:
        print('Querying Gaia... This may take a moment')
        comp_idx = query_gaia.find_comp_bycol2(ref, target_no)
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
        plt.savefig('PSF_Model')
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
