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
from astropy.coordinates import SkyCoord, match_coordinates_sky, ICRS, FK5
from astropy import wcs
from astropy.wcs.utils import pixel_to_skycoord
from astropy.stats import SigmaClip
from astroquery.vizier import Vizier
import pandas as pd
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def exponential_func(x, a, b, c):
        return a * np.exp(b * x) + c
    
def inverse_exponential_func(y, a, b, c):
    if a == 0 or b == 0 or (y - c) / a <= 0:
        return np.nan 
    return (1 / b) * np.log((y - c) / a)

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

def read_sex(sex_file, im):
    """
    A file to read in SExtractor .cat files into a pandas dataframe
    
    params:
    sex_file - The SExtractor catolog path 
    
    return:
    cat_df - The catolog in Pandas Dataframe format
    """
    cat_df = ascii.read(sex_file, format='sextractor')
    im_header = fits.getheader(im)
    im_wcs = wcs.WCS(im_header) 
    
    # Clean spurious detections
    cat_df = cat_df[cat_df['ALPHA_J2000'] != 0.0]
    cat_df = cat_df[cat_df['DELTA_J2000'] != 0.0]
    cat_df = cat_df[cat_df['MAG_APER'] < 0.0]

    eq_coords = pixel_to_skycoord(cat_df['X_IMAGE'] - 1 ,cat_df['Y_IMAGE'] - 1, im_wcs)
    
    # Rename columns to be consistent 
    names = ('NUMBER', 'ALPHA_J2000', 'DELTA_J2000')
    new_names = ('ID', 'RA', 'DEC')
    cat_df.rename_columns(names, new_names)

    cat_df['RA'] = eq_coords.ra.deg
    cat_df['DEC'] = eq_coords.dec.deg
    
    # Convert to df
    cat_df = cat_df.to_pandas()
    
    return cat_df

def cross_match(cat_df, gaia_df, dist):
    """
    A function to cross match two catalogs together using sky coordinates, using astropy
    
    Input:
    
    cat_df - The input photometric catalogue (as a DataFrame) from your telescope (i.e H50)
    gaia_cat - Your reference catalogue (as a DataFrame), default to GAIA DR3
    dist - The max distance to consider a source a match (in arcseconds)
    """
    coord_cat = SkyCoord(ra=cat_df['RA'], dec=cat_df['DEC'], unit=('deg', 'deg'), frame='fk5')
    coord_j2016 = coord_cat.transform_to(FK5(equinox='J2016.0'))
    coord_icrs = coord_cat.transform_to(ICRS)
    coord_gaia = SkyCoord(ra=gaia_df['RA_ICRS'], dec=gaia_df['DE_ICRS'], unit=('deg', 'deg'), frame='icrs')
    
    # Use astropy to match coordinates
    idx, d2d, d3d = match_coordinates_sky(coord_icrs, coord_gaia)
    idx = idx[d2d.arcsec <= float(dist)]
    cat_df['DISTANCE'] = d2d.arcsec
    
    # Clean the two catologs to only sources that match within [dist] arcsec
    cat_df['RA'] = coord_icrs.ra.deg
    cat_df['DEC'] = coord_icrs.dec.deg
    cat_match = cat_df[cat_df['DISTANCE'] <= float(dist)]
    gaia_match = gaia_df.iloc[idx]
    
    return cat_match.reset_index(), gaia_match.reset_index()

def plot_zp_calibration(im,tele_cat, gaia_cat, fltr, phot_dir):
    """
    A function to plot the zeropoint calibration
    
    Input:
    
    tele_cat - The input photometric catalogue from your telescope (i.e H50)
    gaia_cat - Your reference catalogue, default to GAIA DR3
    fltr - The filter of observation, used to find the appropriate values in refernce catalogue
    """
    # Calculate the zeropoint as the difference between the GAIA magnitude 
    # and the instrumental H50 magnitude (MAG_APER), using sigma clipping to clean outliers.
    name = os.path.basename(im)
    mag_dif = gaia_cat[fltr+'mag']-tele_cat['MAG_APER']
    sig = SigmaClip(sigma=2.5, maxiters=10, cenfunc = 'median', stdfunc = 'mad_std')
    clipped_data = sig(mag_dif, masked = True)
    mask = np.array(~clipped_data.mask)  # invert mask (True = keep)
    clipped_mag = clipped_data[mask].data
    clipped_tele = tele_cat[mask]
    clipped_gaia = gaia_cat[mask]
    
    # Calculate statistics for sigma clipped data
    zp_mean, zp_std = np.mean(clipped_data), np.std(clipped_data)
    N = len(clipped_data[clipped_data.mask == False])

    # Figure with histogram on the side
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4)

    # Create the main scatter plot
    ax_scatter = fig.add_subplot(gs[1:4, 0:3])
    ax_scatter.scatter(gaia_cat[fltr+'mag'], mag_dif, color = 'cornflowerblue',\
                       alpha = 0.2, s = 5, label = 'Outliers')
    ax_scatter.errorbar(x = clipped_gaia[str(fltr+'mag')].values, y = clipped_mag, yerr = clipped_tele['MAGERR_APER'].values,\
                        fmt = '.', color= 'blue', alpha = 0.5, label = 'Good Stars (N={:})'.format(N))
    ax_scatter.hlines(zp_mean, gaia_cat[fltr+'mag'].min() - 0.2, gaia_cat[fltr+'mag'].max() + 0.2,\
                      linestyle = 'dashed', color = 'r')
    ax_scatter.grid(alpha = 0.1)
    ax_scatter.set_xlabel("GAIA {:} [mag]".format(fltr), fontsize = 14)
    ax_scatter.set_ylabel('Zeropoint [mag]', fontsize = 14)
    ax_scatter.legend(loc = 'upper right')
    ax_scatter.set_ylim(np.min(mag_dif) - 0.05, np.max(mag_dif) + 0.05)
    ax_scatter.set_xlim(gaia_cat[fltr+'mag'].min() - 0.2, gaia_cat[fltr+'mag'].max() + 0.2)
    ax_scatter.set_title('Zeropoint Calibration - {:}'.format(name), fontweight = 'bold')

    # Display mean values as text on the plot
    ax_scatter.text(0.05, 0.95, 'Mean Zeropoint [mag]: {:.3f}'.format(zp_mean),\
                        transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')
    ax_scatter.text(0.05, 0.90, '$\sigma_{Zp}$ [mag]: '+str(np.round(zp_std/np.sqrt(N),3)),\
                        transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')

    # Create the marginal histogram axes
    ax_histy = fig.add_subplot(gs[1:4, 3], sharey=ax_scatter)

    # Plot the y-axis histogram
    ax_histy.hist(mag_dif, bins=10, orientation='horizontal', color='cornflowerblue', alpha=0.5)
    ax_histy.set_xlabel('Count')
    ax_histy.tick_params(axis='y', labelleft=False)  # Hide y-axis labels for marginal plot
    ax_histy.grid(alpha = 0.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(phot_dir,'Zeropoint_{:}.png'.format(name.replace('.fits',''))))

def get_zp(im, cat_match, gaia_match, fltr, phot_dir, plot = False):
    """
    A function to obtain the magnitude zeropoint, defined as:
    Zp = < m_calibrated - m_instrumental >
    
    Uses sigma clipping to improve estimation, by cutting out bright and faint stars.
    
    """
    mag_dif = gaia_match[fltr+'mag']-cat_match['MAG_APER']
    
    # Perform sigma clipping with sigma = 3 and 10 iterations 
    sig = SigmaClip(sigma=2.5, maxiters=10, cenfunc = 'median', stdfunc = 'mad_std')
    clipped_data = sig(mag_dif, masked = False)
    
    # Calculate statistics for sigma clipped data
    zp_mean, zp_std = np.mean(clipped_data), np.std(clipped_data)
    N = len(clipped_data)

    if plot == True:
        plot_zp_calibration(im, cat_match, gaia_match, fltr, phot_dir)
    
    return zp_mean, zp_std, N

def plot_astrometric_residuals(im, dec_res, ra_res, phot_dir):
    """
    A function to calculate and plot the astrometric residuals
    
    Input:
    
    dec_res - The residuals in declination (arcsec)
    ra_res - The residuals in right ascension (arcsec)
    """
    # Compute median (centroid) and standard deviations for RA and Dec residuals
    name = os.path.basename(im)
    centroid_x = np.median(dec_res)
    dec_std = np.std(dec_res)
    centroid_y = np.median(ra_res)
    ra_std = np.std(ra_res)

    # Create figure, with histograms on the sides
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4)

    # Create the main scatter plot axes
    ax_scatter = fig.add_subplot(gs[1:4, 0:3])
    ax_scatter.scatter(dec_res, ra_res, s = 4, label = 'Stars (N= {:})'.format(len(dec_res)))
    ax_scatter.set_xlabel('Dec Residuals ["]', fontsize = 14)
    ax_scatter.set_ylabel('RA Residuals ["]', fontsize = 14)
    ax_scatter.scatter(centroid_x, centroid_y, marker = 'x', c = 'r', label = 'Centroid')
    ax_scatter.hlines(0, np.min(dec_res) - 0.1, np.max(dec_res) + 0.1,\
                      linestyle = 'dashed', color = 'k', alpha = 0.5)
    ax_scatter.vlines(0, np.min(ra_res) - 0.1, np.max(ra_res) + 0.1, \
                      linestyle = 'dashed', color = 'k', alpha = 0.5)
    ax_scatter.set_title('Astrometric Residuals - {:}'.format(name), fontweight = 'bold')
    ax_scatter.grid(alpha = 0.1)
    ax_scatter.set_xlim(np.min(dec_res) - 0.1, np.max(dec_res) + 0.1)
    ax_scatter.set_ylim(np.min(ra_res) - 0.1, np.max(ra_res) + 0.1)
    ax_scatter.legend(loc = 'upper right')
    
    # Display mean values as text on the plot
    ax_scatter.text(0.05, 0.25, 'Median Dec Offset ["]: {:.3f}'.format(centroid_x),\
                    transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')
    ax_scatter.text(0.05, 0.20, '$St. Dev_{Dec}$ ["]: '+str(np.round(dec_std,3)),\
                    transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')
    ax_scatter.text(0.05, 0.15, 'Median RA Offset ["]: {:.3f}'.format(centroid_y), \
                    transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')
    ax_scatter.text(0.05, 0.10, '$St. Dev_{RA}$ ["]: '+str(np.round(ra_std,3)),\
                    transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')

    # Create the marginal histogram axes
    ax_histx = fig.add_subplot(gs[0, 0:3], sharex=ax_scatter)
    ax_histy = fig.add_subplot(gs[1:4, 3], sharey=ax_scatter)

    # Plot the Dec (x-axis) histogram
    ax_histx.hist(dec_res, bins=10, color='cornflowerblue', alpha=0.7)
    ax_histx.set_ylabel('Count')
    ax_histx.tick_params(axis='x', labelbottom=False)  # Hide x-axis labels for marginal plot
    ax_histx.grid(alpha = 0.1)

    # Plot the RA (y-axis) histogram
    ax_histy.hist(ra_res, bins=10, orientation='horizontal', color='cornflowerblue', alpha=0.7)
    ax_histy.set_xlabel('Count')
    ax_histy.tick_params(axis='y', labelleft=False)  # Hide y-axis labels for marginal plot
    ax_histy.grid(alpha = 0.1)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.savefig(os.path.join(phot_dir,'Astrometric_Residuals_{:}.png'.format(name.replace('.fits',''))))

def get_astrometric_error(im,tele_match, gaia_match, phot_dir, plot = False):
    """
    A function calculate the astrometric errors for an image
    
    """
    # Calculate the residuals and convert from degrees to arcseconds
    ra_res = (gaia_match['RA_ICRS'] - tele_match['RA'])*3600
    dec_res = (gaia_match['DE_ICRS'] - tele_match['DEC'])*3600

    # Compute median (centroid) and standard deviations for RA and Dec residuals
    centroid_dec = np.median(dec_res)
    dec_std = np.std(dec_res)
    centroid_ra = np.median(ra_res)
    ra_std = np.std(ra_res)

    if plot == True:
        plot_astrometric_residuals(im,dec_res, ra_res, phot_dir)

    return centroid_dec, centroid_ra, ra_std, dec_std
    

def plot_maglim(im,tele_cat, param, lim_mag_10, lim_mag_5, fltr, phot_dir):
    """
    Function to plot the error as a function of magnitude and find the
    limiting magnitude at SNR = 10 and 5 (roughly).
    
    Uses the approximation: SNR ~ 1/mag_err
    
    Input:
    
    tele_cat - The calibrated input photometric catalogue from your telescope (i.e H50) 
    lim_mag_10, lim_mag_5 - The limiting magnitudes for SNR = 10 and 5. respectively
    fltr - The filter of observation, used to find the appropriate values in refernce catalogue
    """
    name = os.path.basename(im)
    mag_vec = np.arange(10,21,0.1)
    err_vec = exponential_func(mag_vec, *param)

    # Figure with residual plot on the bottom
    fig = plt.figure(figsize=(8, 6))
    gs = GridSpec(2, 1, height_ratios=[2, 1], hspace=0.1)

    # Create the main scatter plot axes
    ax_scatter = fig.add_subplot(gs[0])
    ax_scatter.scatter(tele_cat['H50_{:}'.format(fltr)], tele_cat['H50_{:}_ERR'.format(fltr)],\
                       label = 'Stars (N = {:})'.format(len(tele_cat)), color = 'blue', alpha = 0.5)
    ax_scatter.plot(mag_vec, err_vec, alpha = 0.5, color = 'k', label = 'Error Model', linestyle = 'dashed')
    ax_scatter.set_ylim(0, 0.105)
    ax_scatter.set_xlim(tele_cat['H50_{:}'.format(fltr)].min() - 0.1, 20)
    ax_scatter.grid(alpha = 0.1)
    ax_scatter.set_ylabel('Error [mag]', fontsize = 14)
    ax_scatter.hlines(0.1, tele_cat['H50_{:}'.format(fltr)].min() - 0.1, 21,\
                      label = 'SNR = 10', color = 'green', linestyle = 'dashed')
    ax_scatter.tick_params(axis='x', which='both',labelbottom=False) 
    ax_scatter.legend(loc = 'upper right')

    # Add the limiting magnitudes as text
    ax_scatter.text(0.05, 0.9, 'Limiting Magnitude (SNR = 10): {:.1f}'.format(lim_mag_10),\
                        transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')
    ax_scatter.text(0.05, 0.80, 'Limiting Magnitude (SNR ~ 5): {:.1f}'.format(lim_mag_5),\
                        transform=ax_scatter.transAxes, fontsize=10, verticalalignment='top')

    # Plot the residuals
    ax_res = fig.add_subplot(gs[1], sharex=ax_scatter)
    ax_res.scatter(tele_cat['H50_{:}'.format(fltr)],\
                   tele_cat['H50_{:}_ERR'.format(fltr)] - exponential_func(tele_cat['H50_{:}'.format(fltr)],*param),\
                   color='blue', alpha=0.5)
    ax_res.set_ylabel('Residuals', fontsize =14)
    ax_res.set_xlabel('r [mag]', fontsize = 14)
    ax_res.grid(alpha = 0.1)
    ax_res.set_ylim(-0.002, 0.002)
    ax_res.hlines(0, 10, 20, linestyle = 'dashed', color = 'k', alpha = 0.5)
    ax_scatter.set_title('Limiting Magnitude - {:}'.format(name), fontweight = 'bold')

    plt.savefig(os.path.join(phot_dir,'LimitingMagnitude_{:}.png'.format(name.replace('.fits',''))))

def calc_maglim(im, tele_cat, fltr, phot_dir, plot = False):
    """
    A function to find the rough magnitude limit for an image at SNR = 10 and SNR = 5.
    """
    # Fit an exponential function to the errors
    param, pcov = curve_fit(exponential_func, tele_cat['H50_{:}'.format(fltr)],\
                            tele_cat['H50_{:}_ERR'.format(fltr)], p0=(0, 0.1, 0.1))
    
    a_fit, b_fit, c_fit = param

    # Calculate limiting magnitude for SNR = 10 and SNR = 5
    Err_SNR_10 = 0.1
    lim_10 = inverse_exponential_func(Err_SNR_10, a_fit, b_fit, c_fit)
    
    Err_SNR_5 = 0.2
    lim_5 = inverse_exponential_func(Err_SNR_5, a_fit, b_fit, c_fit)
    
    # If the detailed method fails, let us approximate with closest data points
    if lim_10 == np.nan or lim_5 == np.nan:
        cat_sub = tele_cat[tele_cat['MAGERR_APER'] <= 0.1]
        lim_10 = cat_sub['H50_{:}'.format(fltr)].max()
        
        cat_sub = tele_cat[tele_cat['MAGERR_APER'] <= 0.2]
        lim_5 = cat_sub['H50_{:}'.format(fltr)].max()

    if plot == True:
        plot_maglim(im, tele_cat,param, lim_10, lim_5, fltr, phot_dir)

    try:
        fits.setval(im, keyword='LIM_10', value=lim_10)
        fits.setval(im, keyword='LIM_5', value=lim_5)
    except:
        raise
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

def add_asterr2hdr(im, ra_err, dec_err):
    """
    A convenience function to add zeropoints to a fits image header.
    
    params:
    im - The fits image 
    ra_res - The standard deviation of RA residuals (arcsec)
    dec_err - The standard deviation of Dec residuals (arcsec)
    """
    fits.setval(im, keyword='RA_ER', value=float(ra_err))
    fits.setval(im, keyword='DEC_ER', value=float(dec_err))


def calibrate_phot(im, cat_df, fltr, phot_dir, plot = False):
    """
    Function to do the work of calibrating a photometry catolog to GAIA DR2.
    Then saves a new catolog in .csv format.
    """
    gaia_df = find_gaia(cat_df, phot_dir)
    
    try:
        # Match catologs and get zeropoints and uncertainty
        cat_match, gaia_match = cross_match(cat_df, gaia_df, 1.0)
        zp, zp_std, N = get_zp(im,cat_match, gaia_match, fltr, phot_dir, plot)
        offset_dec, offset_ra, ra_std, dec_std = get_astrometric_error(im,cat_match, gaia_match, phot_dir, plot)

        # Add zeropoint to the .fits headers
        add_zp2hdr(im, zp, zp_std, N)
        add_asterr2hdr(im, ra_std, dec_std)
    
        # Make a new catolog with calibrated mags
        cat_df['H50_'+fltr] = cat_df['MAG_APER'] + zp
        cat_df['H50_'+fltr+'_ERR'] = np.sqrt((cat_df['MAGERR_APER']**2) + (zp_std/np.sqrt(N))**2)
        calc_maglim(im, cat_df, fltr, phot_dir, plot)
        return cat_df
    
    except:
        raise
        print('Unable to calibrate photometry!')
        raise SystemExit(1) 

        
# Can be called on individual catalogues 
if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Calibrate a photometry catalogue to GAIA DR2 synthetic photometry')
    parser.add_argument('path2cat', type=str, help='The path to the catalogue.')
    parser.add_argument('path2im', type=str, help='The path to the image that the catalogue was extracted from.')
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
        # Match catalogue and get zeropoints and uncertainty
        cat_match, gaia_match = cross_match(cat_df, gaia_df, 1.0)
        zp, zp_std, N = get_zp(im,cat_match, gaia_match, fltr, par_dir, True)
        offset_dec, offset_ra, ra_std, dec_std = get_astrometric_error(im,cat_match, gaia_match, par_dir, True)
    
        # Make a new catalogue with calibrated mags
        cat_df['H50_'+fltr] = cat_df['MAG_APER'] + zp
        cat_df['H50_'+fltr+'_ERR'] = np.sqrt((cat_df['MAGERR_APER']**2) + (zp_std/np.sqrt(N))**2)
        cat_df.to_csv(args.path2cat.replace('.cat','_phot.csv'))
        calc_maglim(im, cat_df, fltr, phot_dir, True)
    except:
        print('Unable to calibrate photometry!')
        raise SystemExit(1) 
    
    
    
    
    
    
    
