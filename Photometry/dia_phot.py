"""
dia_phot.py 

Perform image subtraction photometry using AstroAlign, OIS and PhotUtils packages on H50 images (now in parallel!)

Author: Tom Plunkett

Date: Wed 9 Oct 2024

Organisation: UTAS

"""
# Import necessary package
from prose import Image, Sequence, blocks
import glob
import numpy as np
import argparse
import os
import astroalign as aa
from astropy.io import fits
from astropy.time import Time
from astropy import wcs
from astropy.wcs.utils import pixel_to_skycoord
from astropy.visualization import LinearStretch, ZScaleInterval, ImageNormalize 
import sep
from ois import optimal_system, eval_adpative_kernel
from scipy import ndimage
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
from photutils.centroids import centroid_sources, centroid_com
import pandas as pd
import matplotlib.pyplot as plt
import time
import datetime
import multiprocessing as mp
from multiprocessing import Pool
from calibrate_phot import *

# Let's define any necessary functions

def find_ref(target_dir, norm_ims):
    """ 
    Finds the best reference image (or images) to use
    
    params:
    
    target_dir - The path to the directory with fits files
    
    returns:
    
    best_im - Path to the best image
    best_df - The best night of data in dataframe
    """
    night_dir = os.path.join(target_dir, 'Nightly_Summaries')
    if os.path.isdir(night_dir):
        sum_files = glob.glob(os.path.join(night_dir,'Summary_*.csv'))
        
        for i in range(0, len(sum_files)):
            df = pd.read_csv(sum_files[i])
            df = df.sort_values(['Bkg [adu/s]', 'FWHM [pix]']).reset_index()
            temp = df[df.index == 0]
            
            # Write FWHM/Bkg data to header
            for fl in df.File:
                new_fl = os.path.join(target_dir, os.path.basename(fl))
                fits.setval(new_fl, keyword = 'FWHM',\
                            value = float(df[df.File == fl]['FWHM [pix]'].to_list()[0]), after = 'SEEING')
                fits.setval(new_fl, keyword = 'Bkg',\
                            value = float(df[df.File == fl]['Bkg [adu/s]'].to_list()[0]), after = 'SEEING')
    
            if i == 0:
                best_im = temp
                best_df = df
            if float(temp['FWHM [pix]'][0]) < float(best_im['FWHM [pix]'][0]):
                best_im = temp
                best_df = df 
            
        best_im = os.path.join(target_dir, os.path.basename(best_im['File'][0]))
    else:
        print('Unable to find summary files... Manual mode engaged!')
        fl_choice = str(input('Please enter your desired reference frame (name only): '))
        best_im = os.path.join(target_dir,fl_choice)
        best_df = pd.DataFrame()
    
    return best_im, best_df 

def make_stack(target_dir, best_df, best_im): 
    """ 
    Finds the best reference image (or images) to use
    
    params:
    
    target_dir - The path to the directory with fits files
    
    best_df - The best night summary data in dataframe
    
    returns:
    
    Path to the reference frame file
    """
    data = []
    print(best_df)
    best_im = os.path.join(target_dir, os.path.basename(best_im))
    # If you want a stack, then we will use a simple median
    if len(best_df) > 0 and len(best_df) > 5:
        best_files = best_df['File'][0:5]
    if len(best_df) > 0 and len(best_df) <= 5:
        best_files = best_df['File']
        
    for im in best_files:
        new_fl = os.path.join(target_dir, os.path.basename(im))
        data += [fits.getdata(new_fl.replace('.fits','_reg.fits')).astype(float)]
    
    master_data = np.median(data, axis=0)
    master_header = fits.getheader(best_im.replace('.fits','_reg.fits'))
    master_header['Stack'] = "Median stack of 5 images: {:}".format(best_files.to_list())
    fits.writeto(os.path.join(target_dir, 'ref.fits'), master_data, master_header, overwrite = True)
    
    return os.path.join(target_dir, 'ref.fits')

def get_date(im):
    """
    Get the date of an image from its title (assumes {target_name}_dd_mm_yy_... format)
    
    params:
    
    im - The name of the image
    
    return: 
    date - The string in ddmmyy format
    """
    ar = im.split('_')[1:4]
    date = ar[0] + ar[1] + ar[2]
    return date

def create_subplots(target_dir, norm_ims, diff_ims, positions, ap):
    """
    A function to create a figure to check difference ims against the registered ims
    
    params:
    target_dir - Path to the target directory
    norm_ims, diff_ims - List of paths to the registered and difference images
    pos - The (x,y) tuple of aperture photometry centroid on difference images
    """
    # Define iterators and arrays
    print(ap)
    dates = []
    j = 0
    
    # Get all the dates of observations 
    for im in norm_ims:
        date = get_date(im)
        if date not in dates:
            dates += [date]
    
    # Create the figure in order of date
    fig, ax = plt.subplots(2, len(dates), figsize = (30,10))

    for date in sorted(dates, key=lambda x: datetime.datetime.strptime(x, '%d%m%Y')):
        for i in range(0, len(diff_ims)):
            dt = get_date(diff_ims[i])
            if dt == date:
                norm = fits.getdata(diff_ims[i].replace('_diff', '')).astype(float)
                diff = fits.getdata(diff_ims[i]).astype(float)
                date_obs = fits.getval(diff_ims[i], 'DATE-OBS')
                n1 = ImageNormalize(diff, interval=ZScaleInterval(),stretch=LinearStretch())
                n2 = ImageNormalize(norm, interval=ZScaleInterval(),stretch=LinearStretch())

                ax[0,j].set_title(date_obs)
                ax[0,j].imshow(norm, norm = n2, cmap = 'inferno')
                ax[0,j].set_xlim(250 - 19, 250 + 19)
                ax[0,j].set_ylim(250 - 19, 250 + 19)

                ax[1,j].imshow(diff, norm = n1, cmap = 'inferno')
                ax[1,j].set_xlim(250 - 19, 250 + 19)
                ax[1,j].set_ylim(250 - 19, 250 + 19)
                ax[1,j].add_artist(plt.Circle(positions[0], ap, color='cyan', fill=False, alpha = 1))

                j = j + 1
                break

    plt.tight_layout()
    fig.savefig(os.path.join(target_dir, 'Diff_Check.png'))
           
def onclick(event):
    """
    Event handler for selecting star on image
    """
    global ix, iy 
    if event.dblclick:
        if event.xdata != None and event.ydata != None:
            ix, iy = event.xdata, event.ydata
            
def onclick_multi(event):
    """
    Event handler for selecting comps
    """
    global comp_coords
    if event.dblclick:
        if event.xdata != None and event.ydata != None:
            ix, iy = event.xdata, event.ydata
            comp_coords += [(ix,iy)]
            
def on_esc(event):
    """
    Event handler for selecting comps
    """
    if event.key == 'x' or event.key == 'escape':
        fig.canvas.mpl_disconnect(binding_id)
        plt.close('all')
                       
def register_ims_aa(im_list, ref_im, x,y, verbose):
    """
    A function to register images to the same pixel coordinate grid using AstroAlign
    """
    # Convert to integers for slicing
    x_int = int(x)
    y_int = int(y)
    
    # Load in the reference frame data
    ref_data = fits.getdata(ref_im).astype(float)
    ref_header = fits.getheader(ref_im)
    ref_wcs = wcs.WCS(ref_header)
    ref_mask = ref_data > 30000
    ref_mask = ndimage.binary_dilation(ref_mask, iterations=10)
    masked_ref = np.ma.masked_array(ref_data, ref_mask)
    
    # Iterate through frames in list and register
    for im in im_list:
        if verbose:
            print('Working on frame: {:} \n'.format(im))
        
        # Get the header information and data from current image
        im_header = fits.getheader(im)
        source = fits.getdata(im).astype(float)
        source_wcs = wcs.WCS(im_header)
        source_mask = source > 30000
        source_mask = ndimage.binary_dilation(source_mask, iterations=10)
        masked_source = np.ma.masked_array(source, source_mask)
        
        # Check if WCS is present and remove, as now invalid.
        if source_wcs.has_celestial:
            im_header = im_header[0:(im_header.index('WCSAXES') - 4)]
        
        # If we aren't working on the ref im, use astroalign. Otherwise, just crop to subframe.
        if im != ref_im:
            try:
                registered_image, footprint = aa.register(masked_source, masked_ref,\
                                                          detection_sigma = 10, max_control_points = 200)
                
                registered_image = registered_image[(y_int - 250): (y_int + 250), (x_int - 250): (x_int + 250)]
                new_im = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_reg.fits'))
                fits.writeto(new_im, registered_image, im_header)   
            except:
                print('Failed to register frame: {:}... Moving on! \n'.format(im))
        else:
            source = source[(y_int - 250): (y_int + 250), (x_int - 250): (x_int + 250)]
            new_wcs = ref_wcs[(y_int - 250): (y_int + 250), (x_int - 250): (x_int + 250)]
            ref_header.update(new_wcs.to_header(relax=True))
            new_im = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_reg.fits'))
            fits.writeto(new_im, source, ref_header)
        

def perform_sub(reg_list, ref_im, x,y):
    """
    Function to perform the image subtraction using the Adaptive Bramich algorithm.
    Saves to a new .fits file with _diff.fits at the end. 
    
    params:
    
    reg_list - The list of paths to registered images
    ref_im - The registered ref_im
    
    """
    # Get the reference image data and mask saturated pixels
    ref_data = fits.getdata(ref_im)
    ref_data = np.array(ref_data, dtype = "float")
    ref_mask = ref_data > 30000
    ref_mask = ndimage.binary_dilation(ref_mask, iterations=10)
    ref_bkg = sep.Background(ref_data, mask = ref_mask, bw = 100, bh = 100, fw = 10, fh = 10).back()
    masked_ref = np.ma.masked_array(ref_data - ref_bkg, ref_mask)
    ref_fwhm = fits.getval(ref_im, 'FWHM')
    
    # Iterate through the list
    for im in reg_list:
        if im != ref_im:
            if verbose:
                print('Working on frame: {:} \n'.format(im))
                
            # Get target image data and mask saturated pixels
            im_data = fits.getdata(im).astype('float')
            im_header = fits.getheader(im)
            im_mask = im_data > 30000
            im_mask = ndimage.binary_dilation(im_mask, iterations=10)
            background = sep.Background(im_data, mask = im_mask, bw = 100, bh = 100, fw = 10, fh = 10).back()
            masked_im = np.ma.masked_array(im_data - background, im_mask)
            try:
                im_fwhm = fits.getval(im, 'FWHM')
                fwhm_diff = abs(float(im_fwhm) - float(ref_fwhm))
                kern_size = int(((np.ceil(4*fwhm_diff)//2)*2)+1)
            except:
                kern_size = 13
                
            if kern_size < 7:
                kern_size = 7
            if kern_size > 13:
                kern_size = 13
                
                            
            print(kern_size)
            
            # Subtraction time!
            diff_image, optimal_image, kernel, _ = optimal_system(masked_im, masked_ref,\
                                                               kernelshape = (kern_size,kern_size), bkgdegree = None,\
                                                               method = 'AdaptiveBramich', poly_degree = 2)
            
            # Get the kernel at the object location and find the scale factor
            source_kernel = eval_adpative_kernel(kernel, 250 + (x - int(x)), 250 + (y - int(y)))
            scale_factor = np.sum(source_kernel)
            im_header['x_ob'] = x
            im_header['y_ob'] = y
            im_header['SclFac'] = scale_factor
           
            # Mask bad pixels and save!
            diff_image.data[diff_image.mask] = 0
            fits.writeto(str(im).replace('.fits', '_diff.fits'), diff_image.data, im_header)
            fits.writeto(str(im).replace('.fits', '_kernel.fits'), source_kernel)
            fits.writeto(str(im).replace('.fits', '_bkg.fits'), background)
            
def source_extract(ref_im, target_dir):
    """
    Get a catalog of sources and find zeropoint (using Prose since Sextractor doesn't work anymore)
    
    params: 
    
    ref_im - The path to the reference image
    """
    data = blocks.Get("fluxes", "errors", "exposure", "airmass", 'sky', 'stars_coords', 'wcs', 'filter')
    
    phot = Sequence([
    blocks.detection.PointSourceDetection(),
    blocks.Cutouts(clean = True),                   # making stars cutouts
    blocks.MedianPSF(),                 # building PSF
    blocks.psf.Moffat2D(),              # modeling PSF
    blocks.BalletCentroid(),  
    blocks.PhotutilsAperturePhotometry(apertures = np.array([1.25]), r_in = 2, r_out = 3, scale=True), # aperture photometry
    blocks.Peaks(),
    data,
    ])
    
    # Perform the source extraction
    ref = Image(ref_im)
    phot.run(ref)
    fits.setval(ref_im, keyword = 'FWHM', value = float(ref.fwhm))
    ap = 1.25*float(ref.fwhm)
    
    # Make a dataframe for the source catolog
    sky = np.array(data.sky[0])
    fluxes, fluxes_er = np.array(data.fluxes[0]), np.array(data.errors[0])
    airmass = np.array(data.airmass[0])
    if airmass == None:
        airmass = 0
    mag, mag_er = -2.5*np.log10(fluxes), (2.5/np.log(10))*(fluxes_er/fluxes)
    pos = np.array(data.stars_coords[0])
    WCS = data.wcs[0]
    fltr = str(data.filter[0][0])
    eq_coords = pixel_to_skycoord(pos[:,0], pos[:,1], WCS)
    length = len(fluxes[0])
    d = {'RA': eq_coords.ra, 'DEC': eq_coords.dec, 'MAG_APER': mag[0], 'MAGERR_APER':mag_er[0],\
         'BACKGROUND': [sky]*length, 'X_IMAGE': pos[:,0], 'Y_IMAGE': pos[:,1], 'AIRMASS':[airmass]*length}

    # Calibrate to GAIA
    df = pd.DataFrame(d)
    final_df = calibrate_phot(ref_im, df, fltr, target_dir)
    final_df.to_csv(ref_im.replace('.fits', '_phot.csv'))
    
    return ap
    
def refine_pos(diff_ims):
    """
    Function to refine the position of the object in a simple way
    
    params:
    diff_ims - The list of difference images
    
    return: 
    pos_new - The centroided (x,y) position
    """
    for i in range(0, len(diff_ims)):
        diff = fits.getdata(diff_ims[i]).astype(float)
        if i == 0:
            s = abs(diff)
        if i > 0:
            s = s + abs(diff)
        
    pos_new = centroid_sources(s, 250, 250, box_size = 9, centroid_func=centroid_com)
    pos_new = [(pos_new[0][0], pos_new[1][0])]
    
    print('New centroid at: {:}'.format(pos_new))
    
    return pos_new
     

def target_photometry(target_dir, diff_ims, ref_im, x,y, ap):
    """
    Perform aperture photometry on the difference images for the target object
    
    params:
    diff_ims - List of paths to the difference images
    ref_im - The path to the reference image
    x, y - The coordinates of the target
    
    return: 
    None - Saves the output to a .csv
    """
    # Function to perform the aperture photometry on diff images
    master_df = pd.DataFrame()
    try:
        positions = refine_pos(diff_ims)
    except:
        positions = [(x,y)]
    
    # Create aperture and annulus for bkg estimation
    aperture = CircularAperture(positions, r=ap)
    annulus_aperture = CircularAnnulus(positions, r_in=int(ap*1.2), r_out=int(ap*1.5))
    
    ref_data = fits.getdata(ref_im).astype(float)                                                   
    ref_header = fits.getheader(ref_im)
    date_obs = fits.getval(ref_im, 'DATE-OBS')
    ref_exp = fits.getval(ref_im, 'EXPTIME') 
    
    # Get the zeropoint
    Zp = float(ref_header['MAG_ZP'])
    Zp_er = float(ref_header['ZP_ER'])
    
    ref_table = aperture_photometry(ref_data, aperture)
    ref_bkg = ApertureStats(ref_data, annulus_aperture).median
    ref_flux = ref_table['aperture_sum'] - ref_bkg*aperture.area
    ref_table['Diff. Flux'] = 0
    ref_table['Flux'] = ref_flux
    ref_table['JD'] = float(Time(date_obs, format='isot', scale='utc').jd)
      
    # Iterate through diff images
    for diff in diff_ims:
        scale_factor = float(fits.getval(diff, 'SclFac'))
        data = fits.getdata(diff).astype(float)
        date_obs = fits.getval(diff, 'DATE-OBS')
        try:
            fwhm = fits.getval(diff, 'FWHM')
        except:
            fwhm = np.nan
        exp = fits.getval(diff, 'EXPTIME')
        
        # Do photometry! Suited to 50cm with rn = ... and g = 1.8
        phot_table = aperture_photometry(data, aperture)
        phot_table['Diff. Flux'] = phot_table['aperture_sum']/scale_factor
        phot_table['Flux'] = phot_table['Diff. Flux'] + ref_flux
        phot_table['Flux_Er'] = np.sqrt(phot_table['Flux'] + (aperture.area*(1+1/annulus_aperture.area))*(ref_bkg + (7.10**2) + (1.28**2)/2))
        phot_table['Mag'] = Zp - 2.5*np.log10(phot_table['Flux'])
        phot_table['Mag_Er'] = np.sqrt((Zp_er**2) + ((2.5/np.log(10))*(phot_table['Flux_Er']/phot_table['Flux']))**2)
        phot_table['JD'] = float(Time(date_obs, format='isot', scale='utc').jd)
        phot_table['FWHM'] = fwhm
        phot_table['ExpRatio'] = float(exp)/float(ref_exp)
        phot_table['SclFac'] = scale_factor
        master_df = pd.concat([master_df, phot_table.to_pandas()])
    
    master_df.sort_values('JD').to_csv(os.path.join(target_dir, 'target_dia_phot.csv'))
    
    return positions
                       
if __name__ == '__main__':
    # ---------------------- Part 1 - SETUP ---------------------------------------------------------
     # Hacky shit to stop warnings on mac
    #try:
       # mp.set_start_method('fork', force=True)
   # except RuntimeError:
       # pass
    
    # Set up the parser
    parser = argparse.ArgumentParser(description='Perform DIA analysis on images in folder')
    parser.add_argument('Path', help='The path to folder containing .fits files.')
    parser.add_argument('Visual', type=str, help='Visually choose object? (y/n)')
    parser.add_argument('AutoRef', type=str, help='Find the best reference automatically? (y/n)')
    parser.add_argument('Stack', type=str, help='Stack the best frames? (y/n)')
    parser.add_argument('Redo', type=str, help='Redo subtraction? (y/n)')
    parser.add_argument('Verbose', type=str, help='Want to know whats happening? (y/n)')
    pargs = parser.parse_args()
    
    # Safety
    target_dir = os.path.abspath(pargs.Path)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)
    
    # Define some paths
    dia_dir = os.path.join(target_dir, 'DIA')
    if not os.path.isdir(dia_dir):
        os.mkdir(dia_dir)
    
    # Define some constants and arrays
    ref_auto = str(pargs.AutoRef)
    verbose = str(pargs.Verbose)
    vis = str(pargs.Visual)
    args, args2 = [],[]
    comp_coords = []
    
    if pargs.Redo == 'y':
        # Initiate the pools for multiprocessing
        n_proc = mp.cpu_count()
        print('The number of CPUs is {:}'.format(n_proc))
        p = Pool(int(n_proc/2))
        q = Pool(int(n_proc/2))

        # Remove old files if desired
        old_list = glob.glob(os.path.join(target_dir, '*reg*.fits')) + glob.glob(os.path.join(dia_dir, '*reg*.fits'))
        for old in old_list:
            os.remove(old)

        # Time the run
        start = time.time()

        # Get the paths to the images and choose ref image
        im_list = glob.glob(os.path.join(target_dir, '*_reduced.fits'))

        if ref_auto == 'n' or ref_auto == 'N' or ref_auto == 'No':
            fl_choice = str(input('Please enter your desired reference frame (name only): '))
            ref = os.path.join(target_dir,fl_choice)
        else:
            best_im, best_df = find_ref(target_dir, im_list)

        # Time to pick the coordinates - visual option, else just enter coords.
        ref = Image(best_im)
        if vis == 'y' or vis == 'Y' or vis == 'yes':
            ax = plt.gca()
            fig = plt.gcf()
            ref.show(ax=ax)
            binding_id = fig.canvas.mpl_connect('button_press_event', onclick)
            fig.canvas.mpl_connect('key_press_event', on_esc)
            plt.show()
            x = float(ix)
            y = float(iy)
        else:
            centre_str = input('Input the comma seperated coords. on the reference image (i.e x,y): ')
            x, y = int(centre_str.split(',')[0]), int(centre_str.split(',')[1])

        # ---------------------- Part 2 - Registration and Subtraction ---------------------------------

        # Start image registration in parallel
        print('Performing image registration. Please wait... \n')
        for im in im_list:
            args += [([im], best_im, x, y, 'y')]

        p.starmap(register_ims_aa, args)
        p.close()

        if pargs.Stack == 'y' or pargs.Stack == 'Y' or pargs.Stack == 'Yes':
            ref_im = make_stack(dia_dir, best_df, best_im)
        else:
            ref_im = os.path.join(dia_dir, os.path.basename(best_im).replace('.fits','_reg.fits'))

        ap = source_extract(ref_im, dia_dir)

        # Time to perform subtraction in parallel
        print('Performing image subtraction. Please wait... \n')
        reg_list = glob.glob(os.path.join(dia_dir, '*reg.fits'))

        for im in reg_list:
            args2 += [([im], ref_im, x, y)]

        q.starmap(perform_sub, args2)
        q.terminate()
    
    # ---------------------- Part 3 - Photometry  ---------------------------------------------------------
    else:
        # If they don't want to re-run subtraction, re-choose the ref and proceed
        reg_list = glob.glob(os.path.join(dia_dir, '*reg.fits'))
        
        if ref_auto == 'n' or ref_auto == 'N' or ref_auto == 'No':
            fl_choice = str(input('Please enter the original reference frame chosen (name only): '))
            ref = os.path.join(target_dir,fl_choice)
        else:
            im_list = glob.glob(os.path.join(target_dir, '*_reduced.fits'))
            best_im, best_df = find_ref(target_dir, im_list)
        
        if pargs.Stack == 'y' or pargs.Stack == 'Y' or pargs.Stack == 'Yes':
            ref_im = os.path.join(dia_dir, 'ref.fits')
        else:
            ref_im = os.path.join(dia_dir, os.path.basename(best_im).replace('.fits','_reg.fits'))
        
        x,y = 250,250
        ap = source_extract(ref_im, target_dir)
    
    print('Performing photometry. Please wait... \n')
    diff_list = glob.glob(os.path.join(dia_dir, '*diff.fits'))
    pos = target_photometry(dia_dir, diff_list, ref_im, x, y, ap)
    create_subplots(dia_dir, reg_list, diff_list, pos, ap)
    print('Done!')
