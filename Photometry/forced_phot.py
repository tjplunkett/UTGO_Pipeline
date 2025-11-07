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
from aper_phot import *

# Constants
gain = 1.28 # e-/ADU
rn = 7.10 # e-
D = 0.12 # e/pix/sec

# Define necessary functions
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
parser.add_argument('Path', help='The path to folder containing .fits files.')
parser.add_argument('Aperture', type = float, help = 'The desired aperture for photometry (i.e 5.0)')
parser.add_argument('Visual', type=str, help='Visually choose the target? (y/n)')
parser.add_argument('Verbose', type=str, help='Do you want to know what is going on? (y/n)')
args = parser.parse_args()

# Safety and checks
target_dir = os.path.abspath(args.Path)
if not os.path.isdir(target_dir):
    print("This directory doesn't exist!")
    raise SystemExit(1)
    
#Make an output folder
output_path = os.path.join(target_dir, 'Phot_Outputs')
if not os.path.isdir(output_path):
    os.mkdir(output_path)

r_phot = float(args.Aperture)

# Begin the action!
fm = FitsManager(target_dir, depth=0)

if args.Verbose  == 'Y' or args.Verbose == 'y':
    print(fm)
    
image_list = fm.to_pandas('select path from files')['path'].to_numpy()

# Start! Find object
try:
    ref = Image(fm.all_images[0])
except: 
    ref = Image(image_list[0])

# Choose target position
if args.Visual == 'y':
    # Visual choice
    ax = plt.gca()
    fig = plt.gcf()
    ref.show(ax=ax)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x_obj = float(ix)
    y_obj = float(iy)
    position = [(x_obj, y_obj)]
else:
    x_obj = float(input('Please enter the x coordinate of the target: '))
    y_obj = float(input('Please enter the y coordinate of the target: '))
    position = [(x_obj, y_obj)]

master_phot = pd.DataFrame()

if args.Verbose == 'Y' or args.Verbose == 'y':
    print('Performing image registration and photometry! \n')
    
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
    
    try:
        # Image registration
        transf, (source_list, target_list) = aa.find_transform(ref.data, target_im.data)
        pos = aa.matrix_transform(position, transf.params)

        # Set up the aperture and annulus sizes
        aperture = CircularAperture(pos, r=r_phot)
        annulus_aperture = CircularAnnulus(pos, r_in=3*r_phot, r_out=5*r_phot)

        # Perform aperture photometry on the reference
        phot_table = aperture_photometry(target_im.data, aperture)
        aperstats = ApertureStats(target_im.data, annulus_aperture)
        bkg_mean = aperstats.median
        total_bkg = bkg_mean * aperture.area
        read_noise = ((rn/gain)**2)*aperture.area
        
        # Create a dataframe of the photometry 
        phot_table['JD'] = [jd]
        phot_table['Phot_Aper'] = (phot_table['aperture_sum'] - total_bkg)
        phot_table['SNR'] = phot_table['Phot_Aper']/np.sqrt(phot_table['Phot_Aper']+total_bkg+read_noise)
        phot_table['MAG_APER'] = -2.5*np.log10(phot_table['Phot_Aper'])
        phot_table['MAGERR_APER'] = np.sqrt(((1/phot_table['SNR'])**2))
        phot_table['Airmass'] = [airmass]
        
        if ZP != None and ZP_er != None:
            phot_table['H50_'+str(filt)] = np.round(phot_table['MAG_APER'] + ZP, 3)
            phot_table['H50_'+str(filt)+'_ERR'] = np.round(np.sqrt((phot_table['MAGERR_APER']**2)+(ZP_er**2)),3)
        
        phot_table = phot_table.to_pandas()
        WCS = target_im.wcs
        eq_coords = pixel_to_skycoord(pos[0][0], pos[0][1], WCS)
        phot_table['RA'] = [eq_coords.ra.degree]
        phot_table['DEC'] = [eq_coords.dec.degree]
        master_phot = pd.concat([master_phot, phot_table])
    except:
        print('Unable to perform image registration on {:}... Moving on!'.format(im))

# Save a csv with the photometry 
filename = str(target_object) + '_ap{:}_forcedphot.csv'.format(np.round(r_phot, 2))
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

    # Plot light curve
    plot_lc(master_phot, filt, r_phot, output_path, 'Forced')
