"""
Defringer.py

Date: 15/05/23

Author: Thomas Plunkett

Organisation: UTAS

Purpose:

Scale and subtract a fringe map from a .fits image and then save to a new image.
"""
# Import necessary packages 
import os
import glob
from astropy.io import fits
import numpy as np
import argparse 
import scipy.ndimage as ndimage
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

# Define functions needed
def calc_scalefactor(og_im, fringe_map):
    """
    Calculate the scale factor to multiply the fringe map by
    to match the given image. This is done by minimising the 
    variance in sky background found through sigma clipping
    the difference between the original image and scaled fringe map.
    
    params:
    og_im - The original fringed image
    fringe_map - The fringe map image
    
    returns:
    
    scale_factor_final - The final scale factor (int)
    """
    # Start from 0
    scale_factor = 0
    df_im = og_im - scale_factor*fringe_map
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(df_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    var = np.var(bkg.background)
    
    # The next section is very inefficient and could be improved in future...
    
    # Loop over all scale factors in large steps
    for k in np.arange(0, 10000, 1000):
        df_im = og_im - k*fringe_map
        bkg = Background2D(df_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
        if np.var(bkg.background)<var:
            var = np.var(bkg.background)
            scale_factor = k
    
    # Loop over in smaller steps based upon previous best value
    for j in np.arange(scale_factor - 1000, scale_factor+1000, 100):
        df_im = og_im - j*fringe_map
        bkg_estimator = MedianBackground()
        bkg = Background2D(df_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
        if np.var(bkg.background)<var:
            var = np.var(bkg.background)
            scale_factor = j
    
    # Again, even smaller steps
    for i in np.arange(scale_factor - 100, scale_factor+100, 1):
        df_im = og_im - i*fringe_map
        bkg = Background2D(df_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
        if np.var(bkg.background)<var:
            var = np.var(bkg.background)
            scale_factor_final = i
        
    print('The final scale factor is: ', scale_factor_final)
    return scale_factor_final

def calc_constant(og_im, df_im):
    """
    Calculate a constant to add to the image to match the median count of the de-fringed
    image to the original image (essentially scales back to the initial counts without fringe noise).
    
    params:
    og_im - The original image 
    df_im - The de-fringed image 
    
    return:
    constant - The constant to rescale with
    """
    # Sigma clip and find the background of the original image
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()
    bkg_im = Background2D(og_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    # Sigma clip and find the background of the difference image
    bkg_df = Background2D(df_im, (16, 16), filter_size=(5,5),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
    constant = bkg_im.background_median-bkg_df.background_median
    
    print('The constant is: ',constant)
    return constant

if __name__ == '__main__':
    # Set up the parser
    parser = argparse.ArgumentParser(description='Defringe a set of images')
    parser.add_argument('impath', help='The path to folder containing .fits files.')
    parser.add_argument('mappath', help='The path to the fringe map to be used.')
    args = parser.parse_args()

    target_dir = os.path.abspath(args.impath)
    if not os.path.isdir(target_dir):
        print("This directory doesn't exist!")
        raise SystemExit(1)
       
    # Start the action! 
    map_path = os.path.abspath(args.mappath)

    # Get the fringe image and smooth it using a gaussian filter 
    fringe_im = fits.getdata(map_path)
    fringe_map = fringe_im/np.median(fringe_im)
    fringe_smooth = ndimage.gaussian_filter(fringe_map, sigma=(3,3), order=0, truncate = 2/3)
    
    # Get the fringed images
    im_list = glob.glob(target_dir+'\*.fits')
    for im in im_list:
        print('Working on: ', im)
        image = fits.getdata(im)
        s_f = calc_scalefactor(image, fringe_smooth)
        df_image = image-(s_f*fringe_smooth)
        add = calc_constant(image, df_image)
        final_image = df_image + add
    
        # Write file
        outfile = im.replace('.fits', '_defringed.fits')
        hdu = fits.PrimaryHDU(final_image)
        hdu.header = fits.getheader(im)
        hdu.writeto(outfile, overwrite=True)
    
    print('All done! May the fringe be forever gone...')

                    
