"""
moving_phot.py 

Extract aperture photometry of moving objects from difference images

Author: Tom Plunkett

Date: Friday 7 November 2025

Organisation: UTAS

"""

# Import necessary packages
import glob
import argparse 
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from prose import Image, Sequence, blocks, FitsManager
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
from photutils.centroids import centroid_sources, centroid_com
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from dia_phot import *

def make_skip_stack(target_dir, best_df, best_im): 
    """ 
    Creates a median of best images (to avoid residual signal of moving object)
    
    params:
    
    target_dir - The path to the directory with fits files
    
    best_df - The best night summary data in dataframe
    
    returns:
    
    Path to the reference frame file
    """
    data = []
    print(best_df)
    best_im = os.path.join(target_dir, os.path.basename(best_im))
    
    # If you want a stack, then we will use a simple median on every third image
    if len(best_df) > 0:
        best_files = best_df['File'][::3]
       
    for im in best_files:
        new_fl = os.path.join(target_dir, os.path.basename(im))
        data += [fits.getdata(new_fl.replace('.fits','_reg.fits')).astype(float)]
    
    master_data = np.median(data, axis=0)
    master_header = fits.getheader(best_im.replace('.fits','_reg.fits'))
    master_header['Stack'] = "Median stack of {:} images: {:}".format(len(best_files), best_files.to_list())
    fits.writeto(os.path.join(target_dir, 'ref.fits'), master_data, master_header, overwrite = True)
    
    return os.path.join(target_dir, 'ref.fits')


def register_ims_ast(im_list, ref_im, x_ref,y_ref, verbose):
    """
    A function to register images to the same pixel coordinate grid using AstroAlign
    """
    # Convert to integers for slicing
    x_int = int(x_ref)
    y_int = int(y_ref)
    
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
                raise
                print('Failed to register frame: {:}... Moving on! \n'.format(im))
        else:
            source = source[(y_int - 250): (y_int + 250), (x_int - 250): (x_int + 250)]
            new_wcs = ref_wcs[(y_int - 250): (y_int + 250), (x_int - 250): (x_int + 250)]
            ref_header.update(new_wcs.to_header(relax=True))
            new_im = os.path.join(dia_dir, os.path.basename(im).replace('.fits', '_reg.fits'))
            fits.writeto(new_im, source, ref_header)

def onclick(event):
    """
    Event handler for selecting star on image
    """
    global ix, iy 
    if event.dblclick:
        if event.xdata != None and event.ydata != None:
            ix, iy = event.xdata, event.ydata
            
def on_esc(event):
    """
    Event handler for selecting comps
    """
    if event.key == 'x' or event.key == 'escape':
        fig.canvas.mpl_disconnect(binding_id)
        plt.close('all')

def track_sources_v2(sources_list, initial_position, max_distance=10.0):
    """
    Track a moving object across frames using nearest neighbor matching and velocity-based prediction.

    Parameters:
    - sources_list: List of lists of detected sources (each with `.coords`) in each frame.
    - initial_position: Tuple or list representing the initial (x, y) position of the object.
    - max_distance: Maximum allowed distance to consider a match.

    Returns:
    - tracked_positions: List of (x, y) positions for each frame (or None if undetected and untrackable).
    """
    tracked_positions = [np.array(initial_position)]
    prev_pos = np.array(initial_position)
    velocity = np.array([0.0, 0.0])  # Initial velocity is 0

    for i, sources in enumerate(sources_list[1:], start=1):
        best_match = None
        predicted_pos = prev_pos + velocity  # Predict next position
        min_dist = np.inf

        for s in sources:
            pos = np.array(s.coords)
            dist = np.linalg.norm(pos - predicted_pos)
            if dist < min_dist and dist < max_distance:
                min_dist = dist
                best_match = pos

        if best_match is not None:
            velocity = best_match - prev_pos  # Update velocity estimate
            prev_pos = best_match
            tracked_positions.append(best_match)
        else:
            # No match found: extrapolate position using velocity
            prev_pos = predicted_pos
            tracked_positions.append(None)

    return tracked_positions

def moving_photometry(target_dir, diff_ims, ref_im, x_coords, y_coords, object_name, ap):
    """
    Perform aperture photometry on the difference images for the moving object
    
    params:
    diff_ims - List of paths to the difference images
    ref_im - The path to the reference image
    x, y - The coordinates of the target
    
    return: 
    None - Saves the output to a .csv
    """
    master_df = pd.DataFrame()

    # Reference frame information
    ref_header = fits.getheader(ref_im)
    ref_wcs = WCS(ref_header)
    date_obs = fits.getval(ref_im, 'DATE-OBS')
    ref_exp = fits.getval(ref_im, 'EXPTIME') 
    
    # Get the zeropoint and astrometric errors
    Zp = float(ref_header['MAG_ZP'])
    Zp_er = float(ref_header['ZP_ER'])
    ra_er = float(ref_header['RA_ER'])
    dec_er = float(ref_header['DEC_ER']) 
      
    # Iterate through diff images
    for i in range(0,len(diff_ims)):
        diff = diff_ims[i]
        positions = [(x_coords[i], y_coords[i])]

        if not np.isnan(positions).any():
            # Convert pixel to sky coords
            sky_coord = ref_wcs.pixel_to_world(x_coords[i], y_coords[i])
            ra_deg = sky_coord.ra.deg
            dec_deg = sky_coord.dec.deg
        
            # Define the aperture
            aperture = CircularAperture(positions, r=ap)
            annulus_aperture = CircularAnnulus(positions, r_in=int(ap*2.0), r_out=int(ap*3.0))
            
            # Open the kernel
            kernel = fits.getdata(str(diff).replace('_diff.fits', '_kernel.fits'))
            source_kernel = eval_adpative_kernel(kernel, x_coords[i], y_coords[i])
            scale_factor = np.sum(source_kernel)
            
            # Gather diagnostic information
            date_obs = fits.getval(diff, 'DATE-OBS')
            exp = fits.getval(diff, 'EXPTIME')
            try:
                fwhm = fits.getval(diff, 'FWHM')
            except:
                fwhm = np.nan

            # Do photometry! Suited to 50cm with rn = 7.10 and g = 1.8
            data = fits.getdata(diff).astype(float)
            diff_bkg = ApertureStats(data, annulus_aperture).median
            phot_table = aperture_photometry(data, aperture)
            phot_table['RA'] = ra_deg
            phot_table['RA_Er'] = ra_er/3600 # Convert to degrees
            phot_table['DEC'] = dec_deg
            phot_table['DEC_Er'] = dec_er/3600
            phot_table['Diff. Flux'] = phot_table['aperture_sum']
            phot_table['Flux'] = phot_table['Diff. Flux'] - diff_bkg*aperture.area # As this is moving photometry, F_ref = 0
            phot_table['Flux_Er'] = np.sqrt(phot_table['Flux'] + \
                                            (aperture.area*(1+1/annulus_aperture.area)*(diff_bkg + (7.10**2) + (1.28**2)/2)))
            
            # Convert to magnitude with ZP from reference image
            phot_table['Mag'] = Zp - 2.5*np.log10(phot_table['Flux'])
            phot_table['Mag_Er'] = np.sqrt((Zp_er**2) + ((2.5/np.log(10))*(phot_table['Flux_Er']/phot_table['Flux']))**2)
            
            # Additional Info
            phot_table['JD'] = float(Time(date_obs, format='isot', scale='utc').jd)
            phot_table['ExpRatio'] = float(exp)/float(ref_exp)
            phot_table['SclFac'] = scale_factor
            phot_table['FWHM'] = fwhm
            master_df = pd.concat([master_df, phot_table.to_pandas()])

    master_df.sort_values('JD').to_csv(os.path.join(target_dir, str(object_name+'_Ap'+str(np.round(ap,2))+'_phot.csv')))

    return master_df.sort_values('JD')

if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Perform moving object photometry on difference images')
    parser.add_argument('Path', help='The path to folder containing .fits files.')
    parser.add_argument('Visual', type=str, help='Visually choose object? (y/n)')
    parser.add_argument('AutoRef', type=str, help='Make automatic reference frame? (y/n)')
    parser.add_argument('Redo', type=str, help='Redo subtraction? (y/n)')
    parser.add_argument('Distance', type=int, help='Max distance from expected position to consider a match (e.g 10 pix)')
    parser.add_argument('ID', type=str, help='Target name (to put into phot file name)')
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
    detected_sources = []
    good_ims = []

    # --------------------------- Part 1 - Image Subtraction --------------------------------------
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

        p.starmap(register_ims_ast, args)
        p.close()

        ref_im = make_skip_stack(dia_dir, best_df, best_im)
        ap = source_extract(ref_im, dia_dir, True)

        # Time to perform subtraction in parallel
        print('Performing image subtraction. Please wait... \n')
        reg_list = glob.glob(os.path.join(dia_dir, '*reg.fits'))

        for im in reg_list:
            args2 += [([im], ref_im, x, y, 'y')]

        q.starmap(perform_sub, args2)
        q.terminate()
    
    # ---------------------- Part 3 - Reference Phot  ---------------------------------------------------------
    else:
        # If they don't want to re-run subtraction, re-choose the ref and proceed
        reg_list = glob.glob(os.path.join(dia_dir, '*reg.fits'))
        
        if ref_auto == 'n' or ref_auto == 'N' or ref_auto == 'No':
            fl_choice = str(input('Please enter the original reference frame chosen (name only): '))
            ref = os.path.join(target_dir,fl_choice)
        else:
            im_list = glob.glob(os.path.join(target_dir, '*_reduced.fits'))
            best_im, best_df = find_ref(target_dir, im_list)
        
        ref_im = os.path.join(dia_dir, 'ref.fits')
        ap = source_extract(ref_im, dia_dir, True)

    # --------------------------------- Part 4 - Tracking and Photometry ---------------------------------------

    # Get the time ordered list of files
    fm = FitsManager(dia_dir, depth = 0)
    im_df = fm.files(path = True)
    diff_ims = im_df[im_df.path.str.contains('diff.fits')].path.to_numpy()
        
    # Stage 1 - Do detection and centroiding on difference images
    detect = Sequence([
        blocks.detection.PointSourceDetection(threshold = 5), # stars detection
        blocks.COM(),
    ])

    # Find the object on the first frame
    first = Image(diff_ims[0])
    detect.run(first)
    first.show()
    plt.show()
    source_num = int(input('\n Enter the target number (e.g 1): '))
    initial_position =  first.sources[source_num].coords
    print('\n Initial position is: ', initial_position, '\n')

    # Iterate through difference images and extract the moving objects new position/centroid
    for diff in diff_ims[1:]:
        try:
            im = Image(diff)
            detect.run(im, show_progress = False)
            detected_sources.append(im.sources)
            good_ims.append(diff)
        except:
            print('Issue with file {:}... Moving on!'.format(diff))

    tracked_positions = track_sources_v2(detected_sources, initial_position, pargs.Distance)
    x_coords = [pos[0] if pos is not None else np.nan for pos in tracked_positions]
    y_coords = [pos[1] if pos is not None else np.nan for pos in tracked_positions]

    # Show the trajectory
    plt.close('all')
    first.show()
    plt.scatter(x_coords, y_coords, c = 'r', marker = 'x')
    plt.savefig(os.path.join(dia_dir,'{:}_Track.png'.format(str(pargs.ID))))
    plt.show()

    # Perform the photometry
    data_df = moving_photometry(dia_dir, good_ims, ref_im, x_coords, y_coords, str(pargs.ID), ap)
    data_df = data_df[data_df['Mag_Er'] < 0.15]
    plot_lc(data_df, dia_dir)
    print('Moving object photometry has been extracted! Time to call Bruce Willis')
