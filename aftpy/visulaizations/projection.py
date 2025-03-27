"""
Projection Codes
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from skimage.filters import gaussian
import sunpy.map
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime

def orthographic_array_to_sunpy_map(data, B0_deg, L0_deg, rsun_arcsec=960, date=None):
    """
    Convert an orthographic projected array to a SunPy Map with appropriate metadata.

    Parameters
    ----------
    data : 2D np.ndarray
        Orthographic projected image (square, full-disk).
    B0_deg : float
        Solar B-angle (latitude of disk center) in degrees.
    L0_deg : float
        Carrington longitude of central meridian in degrees.
    rsun_arcsec : float
        Solar radius in arcseconds (default: 960").
    date : str or datetime, optional
        Observation time, e.g. '2025-03-26T00:00:00'. Defaults to now.

    Returns
    -------
    sunpy.map.Map
        A SunPy map object with helioprojective WCS.
    """
    ny, nx = data.shape
    cx = (nx - 1) / 2
    cy = (ny - 1) / 2
    scale_arcsec = 2.0*rsun_arcsec/nx

    # Create FITS-WCS header manually
    header = {
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CRPIX1': cx + 1,
        'CRPIX2': cy + 1,
        'CRVAL1': 0.0,
        'CRVAL2': 0.0,
        'CDELT1': scale_arcsec,
        'CDELT2': scale_arcsec,
        'RSUN_OBS': rsun_arcsec,
        'DSUN_OBS': 1.496e13,  # cm
        'HGLT_OBS': B0_deg,
        'HGLN_OBS': 0.0,  # Observer longitude in Stonyhurst; 0 for Earth
        'CRLN_OBS': L0_deg,  # Carrington longitude of disk center
        'DATE-OBS': date,
        'WAVELNTH': 6173,
        'CONTENT': 'Orthographic Projection of Carrington Map',
        'TELESCOP': 'SIMULATED',
        'INSTRUME': 'CUSTOM',
    }

    return sunpy.map.Map((data, header))

def orthographic_projection_carrington(data, latitudes, longitudes, B0_deg, L0_deg, output_size=512):
    """
    Projects a full-sun Carrington map onto an orthographic disk at a given B0 and L0.

    Parameters
    ----------
    data : 2D np.ndarray
        Input Carrington map of shape (n_lat, n_lon), where latitudes go from -90 to +90 degrees.
    latitudes : 1D np.ndarray
        Array of heliographic latitudes corresponding to the data's rows (in degrees).
    longitudes : 1D np.ndarray
        Array of Carrington longitudes corresponding to the data's columns (in degrees).
    B0_deg : float
        Heliographic latitude of the disk center (degrees).
    L0_deg : float
        Carrington longitude of the central meridian (degrees).
    output_size : int, optional
        Size (pixels) of the output square image. Default is 512.

    Returns
    -------
    projected_image : 2D np.ndarray
        Orthographic projection of the input map, with shape (output_size, output_size).
    """
    # Convert degrees to radians
    B0 = np.deg2rad(B0_deg)
    L0 = np.deg2rad(L0_deg)

    # Set up output image grid in normalized coordinates [-1, 1]
    x = np.linspace(-1, 1, output_size)
    y = np.linspace(-1, 1, output_size)
    X, Y = np.meshgrid(x, y)

    # Mask out points outside the unit circle
    R2 = X**2 + Y**2
    visible = R2 <= 1.0

    # Prepare arrays for output coordinates
    Z = np.zeros_like(X)
    Z[visible] = np.sqrt(1.0 - R2[visible])  # Front hemisphere

    # Convert to heliographic coordinates (inverse of projection)
    x = X[visible]
    y = Y[visible]
    z = Z[visible]

    theta = np.arcsin(y * np.cos(B0) + z * np.sin(B0))
    phi = np.arctan2(x, z * np.cos(B0) - y * np.sin(B0)) + L0

    # Convert to degrees and wrap longitudes
    theta_deg = np.rad2deg(theta)
    phi_deg = np.rad2deg(phi) % 360

    # Interpolate from original data
    lat_idx = np.interp(theta_deg, latitudes, np.arange(len(latitudes)))
    lon_idx = np.interp(phi_deg, longitudes, np.arange(len(longitudes)))
    coords = np.vstack([lat_idx, lon_idx])

    projected_image = np.zeros_like(X)
    projected_image[visible] = map_coordinates(data, coords, order=5, mode='wrap')

    return projected_image
