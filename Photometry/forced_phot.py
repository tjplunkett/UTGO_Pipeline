"""
forced_phot.py

Author: Thomas Plunkett

Purpose:

Perform forced aperture photometry on an image, for the case when the source cannot be detected by 
traditional photometry means (i.e SExtractor or Prose)

"""
# Import necessary packages 
from prose import Image, FitsManager, Sequence, blocks
import astroalign as aa
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.wcs.utils import pixel_to_skycoord
from astropy.time import Time
import subprocess as sub

# Constants
gain = 1.28 # e-/ADU
rn = 7.10 # e-
D = 0.12 # e/pix/sec
ki, ki_er = 0.18, 0.02
kr, kr_er = 0.20, 0.02
kg, kg_er = 0.32, 0.01
kV, kV_er = 0.24, 0.01
kB, kB_er = 0.43, 0.01

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

def onclick(event):
    """
    Event handler for selecting star on image
    """
    global ix, iy 
    if event.dblclick:
        if event.xdata != None and event.ydata != None:
            ix, iy = event.xdata, event.ydata   
            fig.canvas.mpl_disconnect(cid)
            plt.close()
    
# Set up the parser
parser = argparse.ArgumentParser(description='Perform forced aperture photometry on an image')
parser.add_argument('path', help='The path to folder containing .fits files.')
parser.add_argument('depth', type=int, help='How many sub folders to search, i.e if 0 will only access given directory.')
parser.add_argument('astr', type=str, help='Has the astrometry been done? (y/n)')
parser.add_argument('vis', type=str, help='Visually choose the target? (y/n)')
args = parser.parse_args()

# Safety and checks
target_dir = os.path.abspath(args.path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)
    
#Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)

# Begin the action!
fm = FitsManager(target_dir, depth=args.depth)
print(fm)
image_list = fm.to_pandas('select path from files')['path'].to_numpy()

# Start! Find object
try:
    ref = Image(fm.all_images[0])
except: 
    ref = Image(image_list[0])

# Choose target position
if args.vis == 'y':
    # Visual choice, for when not using Keck machine
    ax = plt.gca()
    fig = plt.gcf()
    ref.show(ax=ax)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x_obj = float(ix)
    y_obj = float(iy)
    position = [(x_obj, y_obj)]
    r_phot = float(input('Please enter the aperture radius desired: '))
else:
    x_obj = float(input('Please enter the x coordinate of the target: '))
    y_obj = float(input('Please enter the y coordinate of the target: '))
    r_phot = float(input('Please enter the aperture radius desired: '))
    position = [(x_obj, y_obj)]

master_phot = pd.DataFrame()

for im in image_list:
    # Get image and important quantities
    target_im = Image(im)
    exp_time = float(target_im.header['EXPTIME'])
    
    try:
        airmass = float(target_im.header['AIRMASS'])
    except:
        airmass = 0
    
    target_object = target_im.header['OBJECT']
    jd = Time(target_im.header['DATE-OBS'], format='isot', scale='utc').jd
    filt = target_im.header['FILTER']
    ZP, ZP_er = get_zeropoint(im, r_phot)
    k, k_er = get_ext_coeffs(filt)
    
    try:
        # Image registration
        transf, (source_list, target_list) = aa.find_transform(ref.data, target_im.data)
        pos = aa.matrix_transform(position, transf.params)

        # Set up the aperture and annulus sizes
        aperture = CircularAperture(pos, r=r_phot)
        annulus_aperture = CircularAnnulus(pos, r_in=1.5*r_phot, r_out=2*r_phot)

        # Perform aperture photometry on the reference
        phot_table = aperture_photometry(target_im.data, aperture)
        aperstats = ApertureStats(target_im.data, annulus_aperture)
        bkg_mean = aperstats.median
        total_bkg = bkg_mean * aperture.area

        # Create a dataframe of the photometry 
        phot_table['JD'] = [jd]
        phot_table['Phot_Aper'] = gain*(phot_table['aperture_sum'] - total_bkg)
        phot_table['SNR'] = phot_table['Phot_Aper']/np.sqrt(phot_table['Phot_Aper']+gain*total_bkg)
        phot_table['MAG_APER'] = -2.5*np.log10(phot_table['Phot_Aper'])
        phot_table['MAGERR_APER'] = np.sqrt(((1/phot_table['SNR'])**2)+((k_er*airmass)**2))
        phot_table['Airmass'] = [airmass]
        
        if ZP != None and ZP_er != None:
            phot_table['H50_'+str(filt)] = np.round(phot_table['MAG_APER'] + ZP, 3)
            phot_table['H50_'+str(filt)+'_ER'] = np.round(np.sqrt((phot_table['MAGERR_APER']**2)+(ZP_er**2)),3)
        
        phot_table = phot_table.to_pandas()

        # Get the WCS, if it exists
        if args.astr == 'y':
            WCS = target_im.wcs
            eq_coords = pixel_to_skycoord(pos[0][0], pos[0][1], WCS)
            phot_table['RA'] = [eq_coords.ra.degree]
            phot_table['DEC'] = [eq_coords.dec.degree]

        master_phot = pd.concat([master_phot, phot_table])
    except:
        print('Unable to perform image registration on {:}... Moving on!'.format(im))

# Save a csv with the photometry 
filename = str(target_object) + '_forcedphot.csv'
filename = os.path.join(output_path, filename)

if len(master_phot) != 0:
    master_phot.to_csv(os.path.abspath(filename))
    
    # Image cutout to see aperture placement and size
    ax = plt.gca()
    fig = plt.gcf()
    target_im.show(ax=ax)
    aperture.plot(ax=ax, origin=(0,0), color = 'r')
    ax.set_xlim(pos[0][0] - 51, pos[0][0] + 51)
    ax.set_ylim(pos[0][1] - 51, pos[0][1] + 51)
    fig.savefig(filename.replace('forcedphot.csv','aperture.jpg'))