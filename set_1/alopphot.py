"""
alopphot.py

Functions for finding stars and doing aperture photometry for ALOP Problem Set 1 (2026/2027).

The functions below make use of the photutils package (which should thus be installed in your Python environment):
https://photutils.readthedocs.io/en/stable/index.html

For instruction on installing photutils see:
https://photutils.readthedocs.io/en/stable/getting_started/install.html

The function interfaces are based on the previously used sp.py code.

Version 2026.08.17 - Anthony Brown
"""

__all__ = ["aperphot", "findstars"]

import numpy as np
from photutils.aperture import (
    ApertureStats,
    CircularAnnulus,
    CircularAperture,
    aperture_photometry,
)
from photutils.detection import IRAFStarFinder


def findstars(image, threshold, fwhm):
    """
    Detect the stars in a CCD image.

    Parameters
    ----------
    image : int or float 2D array
        CCD image
    threshold : float
        Minimum threshold for detection. Should be 3-4 sigma above background RMS.
    fwhm : float
        FWHM of the 2D Gaussian convolution kernel. Should be the same as the PSF FWHM.

    Returns
    -------
    sources : Astropy Table
        Table with the list of sources and for each source (id, x, y, flux in ADU, sharpness, roundness)

    >>> from astropy.io import fits
    >>> data = fits.getdata("file.fits")
    >>> sources = findstars(data, threshold, fwhm)
    """
    iraffind = IRAFStarFinder(threshold, fwhm)
    sources = iraffind.find_stars(image)
    sources.rename_column("x_centroid", "x")
    sources.rename_column("y_centroid", "y")
    return sources["id", "x", "y", "flux", "sharpness", "roundness"]


def aperphot(image, xc, yc, apradius, skyradii):
    """
    Perform aperture photometry on the stars in a CCD image.

    Parameters
    ----------
    image : int or float 2D array
        CCD image
    xc : float 1D array
        Vector of x coordinates.
    yc : float 1D array
        Vector of y coordinates.
    apradius : float
        The aperture radius in pixels
    skyradii : float 1D array
        Two element vector giving the inner and outer radii (in pixels) to be used for the sky annulus.

    Returns
    -------
    flux : float 1D array
        List of fluxes in ADU.
    sky : float 1D array
        List of sky background values in ADU

    >>> flux, sky = aperphot(im, x, y, 10.0, [12.0,16.0])
    """
    positions = np.stack([xc, yc]).T.reshape(xc.size, 2)

    # Create the circular apertures around each star position
    aperture = CircularAperture(positions, r=apradius)

    # Create the annulus apertures for sky background estimation
    annulus_aperture = CircularAnnulus(positions, r_in=skyradii[0], r_out=skyradii[1])

    # Calculate statistics for the annulus apertures and use the mean as the sky background per pixel
    aperstats = ApertureStats(image, annulus_aperture)
    bkg_mean = aperstats.mean

    # Calculate the total background within the source aperture
    aperture_area = aperture.area_overlap(image)
    total_bkg = bkg_mean * aperture_area

    # Do the aperture photometry
    phot_table = aperture_photometry(image, aperture)

    # Sky background subtraction
    phot_bkgsub = phot_table["aperture_sum"] - total_bkg

    return phot_bkgsub.value, bkg_mean
