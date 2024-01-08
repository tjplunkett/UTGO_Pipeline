"""
define_telescope.py 

Author: Thomas Plunkett

Purpose:

Defines the header keywords and telescope specifications for the 50cm and 1.3m.
This is necessary for all further scripts.
"""
from prose import Telescope

planewave_dict = dict(
    # Name(s)
    # -------
    name = "PlaneWave 50cm",

    # Keywords
    # --------
    keyword_telescope = "TELESCOP",
    keyword_object = "OBJECT",
    keyword_image_type = "IMAGETYP",
    keyword_light_images = "Light Frame",
    keyword_dark_images = "dark",
    keyword_flat_images = "flat",
    keyword_bias_images = "bias",
    keyword_observation_date = "DATE-OBS",
    keyword_exposure_time = "EXPTIME",
    keyword_filter = "FILTER",
    keyword_airmass = "AIRMASS",
    #keyword_fwhm = "FWHM",
    #keyword_seeing = "SEEING",
    keyword_ra = "OBJCTRA",
    keyword_dec = "OBJCTDEC",
    keyword_jd = "JD",
    keyword_bjd = "BJD",
    #keyword_mjd = "MJD-OBS",
    #keyword_flip = "PIERSIDE",
    keyword_observation_time = "LOCALTIM",

    # Units, formats and scales
    # -------------------------
    ra_unit = "hourangle",
    dec_unit = "deg",
    jd_scale = "utc",
    bjd_scale = "utc",
    #mjd = 0,

    # Specs
    # -----
    trimming = (0, 0), # in pixel along y/x
    read_noise = 7.13, # in ADU
    gain = 1.2, # in e-/ADU
    altitude = 646.0, # in meters
    diameter = 0.508, # in meters
    pixel_scale = 0.796, # in arcseconds
    latlong = [-42.4311, 147.2876],
    saturation = 55000, # in ADU
    hdu = 0
)

cdk20_dict = dict(
    # Name(s)
    # -------
    name = "Planewave CDK20",

    # Keywords
    # --------
    keyword_telescope = "TELESCOP",
    keyword_object = "OBJECT",
    keyword_image_type = "IMAGETYP",
    keyword_light_images = "Light Frame",
    keyword_dark_images = "dark",
    keyword_flat_images = "flat",
    keyword_bias_images = "bias",
    keyword_observation_date = "DATE-OBS",
    keyword_exposure_time = "EXPTIME",
    keyword_filter = "FILTER",
    keyword_airmass = "AIRMASS",
    #keyword_fwhm = "FWHM",
    #keyword_seeing = "SEEING",
    keyword_ra = "OBJCTRA",
    keyword_dec = "OBJCTDEC",
    keyword_jd = "JD",
    keyword_bjd = "BJD",
    #keyword_mjd = "MJD-OBS",
    #keyword_flip = "PIERSIDE",
    keyword_observation_time = "LOCALTIM",

    # Units, formats and scales
    # -------------------------
    ra_unit = "hourangle",
    dec_unit = "deg",
    jd_scale = "utc",
    bjd_scale = "utc",
    #mjd = 0,

    # Specs
    # -----
    trimming = (0, 0), # in pixel along y/x
    read_noise = 7.13, # in ADU
    gain = 1.2, # in e-/ADU
    altitude = 646.0, # in meters
    diameter = 0.508, # in meters
    pixel_scale = 0.796, # in arcseconds
    latlong = [-42.4311, 147.2876],
    saturation = 55000, # in ADU
    hdu = 0
)

big_dict = dict(
    # Name(s)
    # -------
    name = "h127",

    # Keywords
    # --------
    keyword_telescope = "TELESCOP",
    keyword_object = "OBJECT",
    keyword_image_type = "IMAGETYP",
    keyword_light_images = "Light Frame",
    keyword_dark_images = "dark",
    keyword_flat_images = "flat",
    keyword_bias_images = "bias",
    keyword_observation_date = "DATE-OBS",
    keyword_exposure_time = "EXPTIME",
    keyword_filter = "FILTER",
    #keyword_airmass = "AIRMASS",
    #keyword_fwhm = "FWHM",
    #keyword_seeing = "SEEING",
    #keyword_ra = "OBJCTRA",
    #keyword_dec = "OBJCTDEC",
    #keyword_jd = "JD",
    #keyword_bjd = "BJD",
    ##keyword_flip = "PIERSIDE",
    keyword_observation_time = "DATE-OBS",

    # Units, formats and scales
    # -------------------------
    ra_unit = "hourangle",
    dec_unit = "deg",
    jd_scale = "utc",
    bjd_scale = "utc",
    mjd = 0,

    # Specs
    # -----
    trimming = (0, 0), # in pixel along y/x
    read_noise = 10, # in ADU
    gain = 1.3, # in e-/ADU
    altitude = 646.0, # in meters
    diameter = 1.0, # in meters
    pixel_scale = 0.352, # in arcseconds
    latlong = [-42.4311, 147.2876],
    saturation = 65535, # in ADU
    hdu = 0
)


telescope1 = Telescope(planewave_dict)
telescope2 = Telescope(cdk20_dict)
telescope3 = Telescope(big_dict)