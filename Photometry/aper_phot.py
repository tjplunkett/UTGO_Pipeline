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

# Constants
gain = 1.28 # e-/ADU
rn = 7.10 # e-
D = 0.12 # e/pix/sec
ki, ki_er = 0.18, 0.02
kr, kr_er = 0.20, 0.02
kg, kg_er = 0.32, 0.01
kV, kV_er = 0.24, 0.01
kB, kB_er = 0.43, 0.01
ZP_ar = []
ZP_er_ar = []

# Define necessary functions
def get_zeropoint(impath, ap_size):
    """
    A function to obtain the zeropoints for the images. If not already in header,
    must call source extraction process 
    """
    im = Image(impath)
    im_dir = os.path.dirname(impath)
    try:
        zp = float(im.header['MAG_ZP'])
        zp_er = float(im.header['ZP_ER'])
    except:
        try:
            command = 'python pseudosex.py ' + im_dir + ' 0 ' + str(ap_size) + ' y y'
            process = sub.Popen([command], shell=True)
            process.wait()
        except:
            # If there are issues with source extractor, use default Prose but warn
            print('WARNING! Source Extractor issue detected! Defaulting to Prose photometry...')
            command = 'python pseudosex.py ' + im_dir + ' 0 ' + str(ap_size) + ' n y'
            process = sub.Popen([command], shell=True)
            process.wait()
        
        # Read from fits header once calculated
        im = Image(impath)
        zp = float(im.header['MAG_ZP'])
        zp_er = float(im.header['ZP_ER'])
        
    return zp, zp_er

def get_ext_coeffs(filtr):
    """
    Returns the approximate extinction coeff. given a filter for the 50cm
    
    params: 
    filtr - The filter of observation (str)
    
    return:
    zp - The magnitude zeropoint
    zp_er - The error on the zeropoint 
    k - The airmass extinction coeff.
    k_er - The error on the extinction coeff.
    """
    if filtr == "i":
        k, k_er = ki, ki_er
    if filtr == "r":
        k, k_er = kr, kr_er
    if filtr == "g":
        k, k_er = kg, kg_er
    if filtr == 'V':
        k, k_er = kV, kV_er
    if filtr == 'B':
        k, k_er = kB, kB_er
    return  k, k_er

# Set up the parser
parser = argparse.ArgumentParser(description='Perform aperture photometry on images within specified directories and output LC and other plots')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('ID', type=str, help='The target ID.')
parser.add_argument('AutoAperture', type=str, help='Use the automatic aperture? (y/n)')
parser.add_argument('verbose', type=str, help='Want to know whats happening? (y/n)')
args = parser.parse_args()

target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)

# Begin the action!
fm = FitsManager(target_dir, depth=args.depth)
good_ims = []
print(fm)

# --------------------------------------------------------------------------------------
# Stage 1 - Build the PSF model for the images and get median FWHM
PSF = Sequence([
    blocks.detection.PointSourceDetection(), # stars detection
    blocks.Cutouts(clean=False), 
    blocks.detection.LimitStars(min=5), # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
])

try:
    ref = Image(fm.all_images[0])

    # Run this and then show us the stars to pick the target
    PSF.run(ref, show_progress=True)
    ref.show()
    plt.show()
    
except:
    ref = Image(fm.all_images[1])
    # Run this and then show us the stars to pick the target
    PSF.run(ref, show_progress=True)
    ref.show()
    plt.show()

# Choose the target and aperture size
target_no = int(input('Which is the target? '))

# First get the FWHM of each image, find the median
fwhm_data = blocks.Get("fwhm")
get_fwhm = Sequence([
    *PSF[0:-1],                   # apply the same calibration to all images
    blocks.psf.Moffat2D(reference=ref),   # providing a reference improve the PSF optimisation
    fwhm_data,
])

for im in fm.all_images:
    if args.verbose == 'y':
        print('Working on {:}'.format(im))
    try:
        get_fwhm.run(im, show_progress = False)
        good_ims += [im]
    except:
        if args.verbose == 'y':
            print('Failed! Moving on.')
        
fwhm = np.array(fwhm_data.fwhm)
med_fwhm = np.median(fwhm)
print('Median FWHM [pix]: ', np.round(med_fwhm, 1))

# --------------------------------------------------------------------------------------
# Stage 2 - Photometry!
phot_data = blocks.Get("jd_utc", "fluxes", "errors", "airmass", "exposure")
aps = np.arange(0.5, 2, 0.1)
aps = np.append(aps, 3)
    
if args.AutoAperture == 'n' or args.AutoAperture == 'N':
    ap_size = float(input('What aperture size would you like? '))
    ap_index = 0

    fixed_photometry = Sequence([
        *PSF[0:-1],
        blocks.detection.LimitStars(min=12),   # discard images not featuring enough stars
        blocks.Twirl(ref.stars_coords, n=12),  # compute image transformation

        # set stars to the reference ones and apply the inverse
        # transformation previously found to match the ones in the image
        blocks.Set(stars_coords=ref.stars_coords),
        blocks.AffineTransform(data=False, inverse=True),

        blocks.BalletCentroid(),                            # stars centroiding
        blocks.PhotutilsAperturePhotometry(apertures=np.array([ap_size]), r_in = 3.5*med_fwhm, r_out=5.0*med_fwhm, scale = False), # aperture photometry

        # Retrieving data from images in a conveniant way
        phot_data,
    ])

auto_photometry = Sequence([
    *PSF[0:-1],
    blocks.psf.Moffat2D(reference=ref),
    blocks.detection.LimitStars(min=12),   # discard images not featuring enough stars
    blocks.Twirl(ref.stars_coords, n=12),  # compute image transformation

    # set stars to the reference ones and apply the inverse
    # transformation previously found to match the ones in the image
    blocks.Set(stars_coords=ref.stars_coords),
    blocks.AffineTransform(data=False, inverse=True),

    blocks.BalletCentroid(), # stars centroiding
    blocks.PhotutilsAperturePhotometry(apertures=aps, r_in = 3.5, r_out=5.0, scale = med_fwhm), # aperture photometry

    # Retrieving data from images in a conveniant way
    phot_data,
])

filt = ref.header['FILTER']
k, k_er = get_ext_coeffs(filt)

if args.AutoAperture == 'n' or args.AutoAperture == 'N': 
    fixed_photometry.run(good_ims)
else: 
    auto_photometry.run(good_ims)
    
# Calculate inst.mag and convert to pandas df
time = np.array(phot_data.jd_utc)
exp = np.array(phot_data.exposure[0].value)
fluxes = np.array(phot_data.fluxes)
fluxes_er = np.array(phot_data.errors)

# Find the best aperture
if args.AutoAperture == 'y' or args.AutoAperture == 'Y': 
    ap_index = np.argmax(fluxes[0, :, target_no]/fluxes_er[0, :, target_no])
    ap_size = med_fwhm*aps[ap_index]
    
airmass = np.array(phot_data.airmass)
flux = fluxes[:, ap_index, target_no]
flux_er = fluxes_er[:, ap_index, target_no]
mag = -2.5*np.log10(flux/exp)
mag_er = (2.5/np.log(10))*(flux_er/flux)

for im in good_ims:
    ZP, ZP_er = get_zeropoint(im, ap_size)
    ZP_ar += [ZP]
    ZP_er_ar += [ZP_er]

cal_mag = -2.5*np.log10(flux) + np.array(ZP_ar)
cal_er = np.sqrt(mag_er**2 + (np.array(ZP_er_ar))**2)
# --------------------------------------------------------------------------------------
# Stage 3 - Save your results!
d = {'JD': time, 'Flux': flux, 'Flux_Er': flux_er, 'Inst_Mag':mag, 'Inst_Mag_Er':mag_er, 'Airmass':airmass, 'Exposure': exp, 'H50_'+str(filt): np.round(cal_mag,3), 'H50_'+str(filt)+'_ER': np.round(cal_er, 3)}
df = pd.DataFrame(d)

# Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
    
# Save file     
file_name = str(os.path.basename(target_dir))+'_aperphot_'+str(args.ID)+'_ap'+str(np.round(ap_size, 2))+'_t'+str(target_no)+'.csv'
df.to_csv(os.path.join(output_path, file_name))