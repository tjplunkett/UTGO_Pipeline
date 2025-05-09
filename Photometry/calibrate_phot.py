"""
calibrate_phot.py

Author: Thomas Plunkett

Purpose: 

Defines functions to calibrate photometry catologs to GAIA DR2, along with finding zeropoints for images.

"""
# Import necessary packages
import argparse
from astropy.io import ascii, fits
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.stats import SigmaClip
from astroquery.vizier import Vizier
import pandas as pd
import os

def find_gaia(cat_df, phot_dir):
    """
    Function to obtain a dataframe containing GAIA DR2 synthetic photometry.
    If cannot find a .csv file, will download using Vizier and save for use later.
    
    params:
    
    cat_df - The dataframe containing the photometry catalog from H50
    """
    if phot_dir != None:
        if os.path.isfile(os.path.join(phot_dir,'gaia_cat.csv')):
            gaia_df = pd.read_csv(os.path.join(phot_dir,'gaia_cat.csv'))
        else:
            v = Vizier(columns=['**'], catalog='I/360/syntphot', row_limit=-1)
            med_ra, med_dec = np.median(cat_df['RA'][:]), np.median(cat_df['DEC'][:])
            coord = SkyCoord(ra=med_ra, dec=med_dec, unit=(u.degree, u.degree), frame='icrs')
            fov = abs(max(cat_df['DEC']) - min(cat_df['DEC'])) * 0.75
            try:
                gaia_cat = v.query_region(coord, radius=fov * u.deg)
                gaia_df = gaia_cat[0].to_pandas()
                gaia_df.to_csv(os.path.join(phot_dir,'gaia_cat.csv'))
            except:
                print('Unable to query GAIA... Please try again!')
                gaia_df = pd.DataFrame() 
            
    else:
        v = Vizier(columns=['**'], catalog='I/360/syntphot', row_limit=-1)
        med_ra, med_dec = np.median(cat_df['RA'][:]), np.median(cat_df['DEC'][:])
        coord = SkyCoord(ra=med_ra, dec=med_dec, unit=(u.degree, u.degree), frame='icrs')
        fov = abs(max(cat_df['DEC']) - min(cat_df['DEC'])) * 0.75
        try:
            gaia_cat = v.query_region(coord, radius=fov * u.deg)
            gaia_df = gaia_cat[0].to_pandas()
            gaia_df.to_csv(os.path.join(phot_dir,'gaia_cat.csv'))
        except:
            print('Unable to query GAIA... Please try again!')
            gaia_df = pd.DataFrame()  
    
    return gaia_df

def read_sex(sex_file):
    """
    A file to read in SExtractor .cat files into a pandas dataframe
    
    params:
    sex_file - The SExtractor catolog path 
    
    return:
    cat_df - The catolog in Pandas Dataframe format
    """
    cat_df = ascii.read(sex_file, format='sextractor')
    
    # Clean spurious detections
    cat_df = cat_df[cat_df['ALPHA_J2000'] != 0.0]
    cat_df = cat_df[cat_df['DELTA_J2000'] != 0.0]
    cat_df = cat_df[cat_df['MAG_APER'] < 0.0]
    
    # Rename columns to be consistent 
    names = ('NUMBER', 'ALPHA_J2000', 'DELTA_J2000')
    new_names = ('ID', 'RA', 'DEC')
    cat_df.rename_columns(names, new_names)
    
    # Convert to df
    cat_df = cat_df.to_pandas()
    
    return cat_df

def cross_match(cat_df, gaia_df):
    """
    A function to cross match two catalogs together using sky coordinates, using astropy
    """
    coord_cat = SkyCoord(ra=cat_df['RA'], dec=cat_df['DEC'], unit=('deg', 'deg'), frame='icrs')
    coord_gaia = SkyCoord(ra=gaia_df['RA_ICRS'], dec=gaia_df['DE_ICRS'], unit=('deg', 'deg'), frame='icrs')
    
    # Use astropy to match coordinates
    idx, d2d, d3d = match_coordinates_sky(coord_cat, coord_gaia)
    idx = idx[d2d.arcsec <= 0.8]
    cat_df['DISTANCE'] = d2d.arcsec
    
    # Clean the two catologs to only sources that match within 0.8 arcsec 
    cat_match = cat_df[cat_df['DISTANCE'] <= 0.8]
    gaia_match = gaia_df.iloc[idx]
    
    return cat_match.reset_index(), gaia_match.reset_index()

def get_zp(cat_match, gaia_match, fltr):
    """
    A function to obtain the magnitude zeropoint, defined as:
    Zp = < m_calibrated - m_instrumental >
    
    Uses sigma clipping to improve estimation, by cutting out bright and faint stars.
    
    """
    mag_dif = gaia_match[fltr+'mag']-cat_match['MAG_APER']
    
    # Perform sigma clipping with sigma = 2 and 5 iterations 
    sig = SigmaClip(sigma=2, maxiters=5)
    clipped_data = sig(mag_dif, masked = False)
    
    # Calculate statistics for sigma clipped data
    zp_mean, zp_std = np.mean(clipped_data), np.std(clipped_data)
    N = len(clipped_data)
    
    return zp_mean, zp_std, N

def calc_maglim(im, cat_df, zp):
    """
    A function to find the rough magnitude limit for an image at SNR = 10 and SNR = 5.
    """
    cat_df['H50_mag'] = cat_df['MAG_APER'] + zp

    # For SNR = 10:
    cat_sub = cat_df[cat_df['MAGERR_APER'] <= 0.1]
    lim_10 = cat_sub['H50_mag'].max()

    # For SNR = 5:
    cat_sub = cat_df[cat_df['MAGERR_APER'] <= 0.2]
    lim_5 = cat_sub['H50_mag'].max()

    try:
        fits.setval(im, keyword='LIM_10', value=lim_10)
        fits.setval(im, keyword='LIM_5', value=lim_5)
    except:
        print('Unable to find the magnitude limits...')
    

def add_zp2hdr(im, zp_mean, zp_std, N):
    """
    A convenience function to add zeropoints to a fits image header.
    
    params:
    im - The fits image 
    zp_mean - The magnitude zeropoint
    zp_std - The standard deviation of zeropoint
    N - Number of stars used to calculate zeropoint
    """
    fits.setval(im, keyword='MAG_ZP', value=zp_mean)
    fits.setval(im, keyword='ZP_ER', value=float(zp_std/np.sqrt(N)))


def calibrate_phot(im, cat_df, fltr, phot_dir):
    """
    Function to do the work of calibrating a photometry catolog to GAIA DR2.
    Then saves a new catolog in .csv format.
    """
    gaia_df = find_gaia(cat_df, phot_dir)
    
    try:
        # Match catologs and get zeropoints and uncertainty
        cat_match, gaia_match = cross_match(cat_df, gaia_df)
        zp, zp_std, N = get_zp(cat_match, gaia_match, fltr)

        # Add zeropoint to the .fits headers
        add_zp2hdr(im, zp, zp_std, N)
        calc_maglim(im, cat_df, zp)
    
        # Make a new catolog with calibrated mags
        cat_df['H50_'+fltr] = cat_df['MAG_APER'] + zp
        cat_df['H50_'+fltr+'_ERR'] = np.sqrt((cat_df['MAGERR_APER']**2) + (zp_std/np.sqrt(N))**2)
        return cat_df
    
    except:
        print('Unable to calibrate photometry!')
        raise SystemExit(1) 

        
# Can be called on individual catalogues 
if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Calibrate a photometry catalogue to GAIA DR2 synthetic photometry')
    parser.add_argument('path2cat', type=str, help='The path to the catalogue.')
    parser.add_argument('filter', type=str, help="The bandpass/filter of observation (i.e for SDSS r', type: r")
    parser.add_argument('sex', type=str, help = 'Is the file a Sextractor .cat file? (y/n)')
    args = parser.parse_args()
    
    # Safety check
    if not os.path.isfile(args.path2cat):
        print('Not a valid file! Check path and try again.')
        raise SystemExit(1)
    
    # Get the parent directory of file
    par_dir = os.path.dirname(args.path2cat)
    
    # Read in catalogue, either a .cat from SExtractor or .csv from PhotUtils
    if args.sex == 'y' or args.sex == 'Y':
        cat_df = read_sex(args.path2cat)
    else:
        cat_df = pd.read_csv(args.path2cat)
    
    # Query gaia and save result for future use
    gaia_df = find_gaia(cat_df, par_dir)
    
    try:
        # Match catologs and get zeropoints and uncertainty
        cat_match, gaia_match = cross_match(cat_df, gaia_df)
        zp, zp_std, N = get_zp(cat_match, gaia_match, args.filter)
        
        print('The zeropoint was determined to be: {:} +/- {:}, using {:} stars'.format(zp, zp_std/np.sqrt(N), N))
    
        # Make a new catolog with calibrated mags
        cat_df['H50_'+args.filter] = cat_df['MAG_APER'] + zp
        cat_df['H50_'+args.filter+'_ERR'] = np.sqrt((cat_df['MAGERR_APER']**2) + (zp_std/np.sqrt(N))**2)
        cat_df.to_csv(args.path2cat.replace('.cat', '_phot.csv'))
    except:
        print('Unable to calibrate photometry!')
        raise SystemExit(1) 
    
    
    
    
    
    
    
