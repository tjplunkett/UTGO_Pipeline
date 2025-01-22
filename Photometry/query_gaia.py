"""
query_gaia.py

Author: Thomas Plunkett

Purpose:

Defines functions to query gaia and find stars of similiar colour and magnitude to a target,
which can be used in differential photometry. Will fail without internet, as
GAIA needs to be queried. 
"""

# Import necessary packages
from prose import Image, Sequence, blocks, FitsManager
from prose.blocks import catalogs

def find_comp_bycol(path2im, source_no):
    """
    A function to find Gaia stars that are closest in colour to the target. This one
    works from calibrated images in a folder, for use outside of diff phot script.
    
    params: 
    path2im - The path to the directory containing calibrated images
    source_no - The label for the source
    
    return:
    comp_index - The index of the stars in the frame to use
    """
    col_lim = 0.2
    
    # Parse the fits files
    fm = FitsManager(path2im, depth=2)
    
    # Use the first image as the reference
    ref = Image(fm.all_images[0])
    
    # Run source extraction
    calibration = Sequence([
        blocks.SEDetection(threshold=3.0), # stars detection
        blocks.Cutouts(),                   # making stars cutouts
        blocks.MedianPSF(),                 # building PSF
        blocks.psf.Moffat2D(),              # modeling PSF
    ])
    
    calibration.run(ref, show_progress=True)
    
    # Query Gaia 
    im = catalogs.GaiaCatalog(mode='crossmatch')(ref)
    im.show(stars=True)
    gaia_raw = im.catalogs['gaia'][['index','bp_rp','phot_bp_mean_mag','x','y', 'phot_variable_flag']]
    gaia_table = gaia_raw[gaia_raw['phot_variable_flag'] != 'VARIABLE']
    print('Flagged {:} stars as variable!'.format(len(gaia_raw) - len(gaia_table)))
    
    # Extract a table for the source, set variables
    source_table = gaia_table.iloc[[source_no]]
    bp_rp = source_table['bp_rp'][source_no]
    bp = source_table['phot_bp_mean_mag'][source_no]
    
    
    # Create absolute residuals for comparison
    gaia_table['col_res'] = abs(gaia_table['bp_rp']-bp_rp)
    gaia_table['mag_res'] = abs(gaia_table['phot_bp_mean_mag'] - bp)
    
    # Make a subset containing only stars with colours within colour limit (i.e starts at 0.2). If the length
    # of this is less than 5 stars, then expand the search by +0.05 mag.
    close_col = gaia_table[gaia_table['col_res']<=col_lim].sort_values('col_res')[1:]
    while len(close_col['mag_res']) < 5:
        print(col_lim)
        col_lim += 0.05
        close_col = gaia_table[gaia_table['col_res']<=col_lim].sort_values('col_res')[1:]
        
    # Take the 5 closest colour and mag stars 
    comp_stars = close_col.sort_values('mag_res')[0:5]
    comp_index = comp_stars.index.to_list()
    
    return comp_index

def find_comp_bycol2(ref, source_no, thresh):
    """
    A function to find Gaia stars that are closest in colour to the target. Works from
    reference image object, to be used within the diff phot script.
    
    params: 
    ref - The reference image object from prose
    source_no - The label for the source
    
    return:
    comp_index - The index of the stars in the frame to use
    
    """
    col_lim = 0.2
    
    # Run source extraction
    calibration = Sequence([
        blocks.detection.PointSourceDetection(threshold = thresh),  # stars detection
        blocks.Cutouts(clean = False), # making stars cutouts
    ])
    
    calibration.run(ref, show_progress=True)
    
    # Query Gaia 
    im = catalogs.GaiaCatalog(mode='crossmatch')(ref)
    im.show(stars=True)
    gaia_raw = im.catalogs['gaia'][['index','bp_rp','phot_bp_mean_mag','x','y', 'phot_variable_flag']]
    x_min,x_max = min(gaia_raw['x']), max(gaia_raw['x'])
    y_min,y_max = min(gaia_raw['y']), max(gaia_raw['y'])
    gaia_table = gaia_raw[gaia_raw['phot_variable_flag'] != 'VARIABLE']
    print('Flagged {:} stars as variable!'.format(len(gaia_raw) - len(gaia_table)))

    # Measure the source colour
    source_table = gaia_table.iloc[[source_no]]
    bp_rp = source_table['bp_rp'][source_no]
    bp = source_table['phot_bp_mean_mag'][source_no]
    x_source = source_table['x'][source_no]
    y_source = source_table['y'][source_no]

    # Remove stars too close to the edge
    #gaia_table = gaia_table[gaia_table['x'] < (x_source + 600)]
    #gaia_table = gaia_table[gaia_table['x'] > (x_source - 600)]
    #gaia_table = gaia_table[gaia_table['y'] < (y_source + 600)]
    #gaia_table = gaia_table[gaia_table['y'] > (y_source - 600)]

    gaia_table = gaia_table[gaia_table['x'] < (x_max - 250)]
    gaia_table = gaia_table[gaia_table['x'] > (x_min + 250)]
    gaia_table = gaia_table[gaia_table['y'] < (y_max - 250)]
    gaia_table = gaia_table[gaia_table['y'] > (y_min + 250)]
    print('Flagged {:} stars as too close to the edge!'.format(len(gaia_raw) - len(gaia_table)))
    print(gaia_table)
   
    # Create absolute residuals for comparison
    gaia_table['col_res'] = abs(gaia_table['bp_rp']-bp_rp)
    gaia_table['mag_res'] = abs(gaia_table['phot_bp_mean_mag'] - bp)
    
    # Make a subset containing only stars with colours within colour limit (i.e starts at 0.2). If the length
    # of this is less than 5 stars, then expand the search by +0.05 mag.
    close_col = gaia_table[gaia_table['col_res']<=col_lim].sort_values('col_res')[1:]
    while len(close_col['mag_res']) < 5:
        print(col_lim)
        col_lim += 0.05
        close_col = gaia_table[gaia_table['col_res']<=col_lim].sort_values('col_res')[1:]
        
    # Take the 5 closest colour and mag stars 
    comp_stars = close_col.sort_values('mag_res')[0:5]
    comp_index = comp_stars.index.to_list()
    
    return comp_index

    
    
    
    
