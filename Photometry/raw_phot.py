"""
raw_phot.py

Author: Thomas Plunkett

Purpose:

Performs aperture photometry on a series of images, in order to get a light curve. This
will be more prone to error than differential photometry.

"""
# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager, Observation
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os

# Set up the parser
parser = argparse.ArgumentParser(description='Perform aperture photometry on images within specified directories and output LC and other plots')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('ID', type=str, help='The target ID.')
parser.add_argument('AutoAperture', type=str, help='Use the automatic aperture? (y/n)')
args = parser.parse_args()

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
    #blocks.Calibration(darks=fm.all_darks, bias=fm.all_bias, flats=fm.all_flats),
    blocks.Trim(),
    blocks.SegmentedPeaks(),             # stars detection
    blocks.Cutouts(),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])

# Run this and then show us the stars to pick the target
calibration.run(ref, show_progress=True)
ref.show()
plt.show()

data = blocks.Get("jd_utc", "fluxes", "errors", "airmass", "exposure", "fwhm")

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
    blocks.PhotutilsAperturePhotometry(apertures=np.arange(0.25, 10.25, 0.25), scale=float(ref.fwhm)), # aperture photometry
    blocks.Peaks(),

    # Retrieving data from images in a conveniant way
    data,

    # Stack image
    blocks.AffineTransform(stars=False, data=True),
    blocks.Stack(ref, name="stack"),
])

# Run the photometry on all images
photometry.run(fm.all_images)

fwhm = np.array(data.fwhm)

# Choose the target and aperture size
target_no = int(input('Which is the target?'))

if args.AutoAperture == 'n':
    ap_size = int(input('What aperture size would you like?'))
    ap_index = int(np.round(((ap_size/ref.fwhm)/0.25),0)) - 1 
    print(ap_index)
else:
    ap_size = np.round(1.5*np.median(fwhm), 0)
    ap_index = int(np.round(((ap_size/ref.fwhm)/0.25),0)) - 1
    print(ap_index)

# Calculate inst.mag and convert to pandas df
time = np.array(data.jd_utc)
exp = np.array(data.exposure[0].value)
fluxes = np.array(data.fluxes)
fluxes_er = np.array(data.errors)
airmass = np.array(data.airmass)
flux = fluxes[:, ap_index, target_no]
flux_er = fluxes_er[:, ap_index, target_no]
mag = -2.5*np.log10(flux/exp)
mag_er = (2.5/np.log(10))*(flux_er/flux)

d = {'JD': time, 'Flux': flux, 'Flux_Er': flux_er, 'Inst_Mag':mag, 'Inst_Mag_Er':mag_er, 'Airmass':airmass, 'FWHM': fwhm, 'Exposure': exp}
df = pd.DataFrame(d)

# Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
    
# Save file     
file_name = str(os.path.basename(target_dir))+'_rawflux_'+str(args.ID)+'_ap'+str(ap_size)+'.csv'
df.to_csv(os.path.join(output_path, file_name))