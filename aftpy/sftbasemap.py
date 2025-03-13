"""
SFTMap Module: A Python Library for Handling AFT Magnetic Field Data

This module provides classes and utilities to read, analyze, and visualize
magnetic field data from the Advanced Flux Telescope (AFT). The data is structured
in FITS format and includes metadata describing solar magnetic flux density.

Classes
-------
- AFTSunpyMap : A subclass of SunPy's GenericMap, designed for handling AFT maps.
- SFTMap : A class to process individual AFT map files, compute various solar magnetic
           parameters, and convert them into SunPy-compatible maps.

Features
--------
- Reads AFT map data from FITS files.
- Provides WCS (World Coordinate System) support for coordinate transformations.
- Computes essential solar magnetic parameters such as:
  - Total absolute flux
  - Axial and equatorial dipole moments
  - Polar region magnetic strengths
- Supports SunPy integration for visualization.
- Provides efficient computation with NumPy-based optimizations.

Usage Example
-------------
>>> from sftmap import SFTMap
>>> sft = SFTMap("path/to/data.fits")
>>> print(sft.parameters())

References
----------
* Advanced Flux Telescope (AFT) Mission Paper
* SunPy Library: https://sunpy.org/
"""

import os
from pathlib import Path
from typing import Union

import numpy as np
import astropy.units as u
from astropy.io.fits import Header
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from sunpy.map import GenericMap
from sunpy.sun.constants import radius
from sunpy.util.metadata import MetaDict

from .io.io import read_aftv2
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from .utilities import detect_sftmodel


class AFTSunpyMap(GenericMap):
    """
    A SunPy-based map for the Next Generation Telescope (AFT).

    Attributes
    ----------
    meta : dict
        Metadata dictionary.
    plot_settings : dict
        Dictionary for visualization settings.
    """

    def __init__(self, data: np.ndarray, header: Union[MetaDict, Header, dict], **kwargs):
        header.setdefault("instrume", "AFT+HMI")
        if header["hgln_obs"] != 0: header["hgln_obs"] =0.0
        super().__init__(data, header, **kwargs)

        self.meta["quantity"] = "Magnetic Flux Density [G]"
        self.meta["obsrvtry"] = "AFTv2.0"
        self.plot_settings["cmap"] = "hmimag"
        self.plot_settings["norm"] = colors.Normalize(vmax=1000, vmin=-1000)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if the provided data and header correspond to an AFT image.

        Parameters
        ----------
        data : np.ndarray
            Image data array.
        header : dict or MetaDict
            FITS header or metadata dictionary.

        Returns
        -------
        bool
            True if the data is from the NextGenerationTelescope, else False.
        """
        return header.get("instrume", "").startswith("AFTv2.0")


class SFTMap:
    """
    A class to handle and analyze AFT map data.

    Parameters
    ----------
    data : np.ndarray or str
        AFT map data or file path.
    header : dict, MetaDict, or Header, optional
        Metadata dictionary.
    mask : np.ndarray, optional
        2D mask array.
    magmap : np.ndarray, optional
        Magnetic field map.

    Attributes
    ----------
    filename : Path or None
        Relative file path.
    data : np.ndarray
        AFT map data.
    header : dict, MetaDict, or Header
        Metadata dictionary.
    mask : np.ndarray
        Mask data.
    magmap : np.ndarray
        Magnetic field map.
    """

    def __init__(self, data: Union[np.ndarray, str], header=None, mask=None, magmap=None, **kwargs):
        self.filename = None

        if isinstance(data, str) and header is None:
            _sftdata = read_aftv2(data)
            _root = Path(data).parent.parent.parent
            self.filename = Path(data).relative_to(_root)

            self.data = _sftdata["data"]
            self.header = _sftdata["header"]
            self.mask = _sftdata["mask"]
            self.magmap = _sftdata["magmap"]
        else:
            self.data = data
            self.header = header
            self.mask = mask
            self.magmap = magmap

    @property
    def shape(self) -> tuple:
        """Returns the shape of the data array."""
        return self.data.shape

    @property
    def wcs(self) -> WCS:
        """Returns the World Coordinate System (WCS) from the header."""
        return WCS(self.header)

    @property
    def timestamp(self) -> str:
        """Extracts the observation timestamp from the header."""
        return self.header.get("date_obs", None)

    @property
    def sunpymap(self) -> AFTSunpyMap:
        """Returns a SunPy-compatible map of the data."""
        return AFTSunpyMap(self.data, self.header)

    @staticmethod
    def _validate_array(data, name):
        if data is None or not isinstance(data, np.ndarray):
            raise ValueError(f"{name} must be a valid NumPy array.")

    @property
    def data(self) -> np.ndarray:
        return self._data

    @data.setter
    def data(self, value):
        self._validate_array(value, "Data")
        self._data = value

    @property
    def header(self) -> Union[MetaDict, Header, dict]:
        return self._header

    @header.setter
    def header(self, value):
        if not isinstance(value, (MetaDict, Header, dict)):
            raise ValueError("Header must be a MetaDict, FITS Header, or dictionary.")
        self._header = value

    @property
    def mask(self) -> np.ndarray:
        return self._mask

    @mask.setter
    def mask(self, value):
        if value is not None and not isinstance(value, np.ndarray):
            raise ValueError("Mask must be a 2D NumPy array.")
        self._mask = value

    def pixel_coord(self, radian: bool = True):
        """
        Computes pixel coordinates in latitude, longitude, and area.

        Parameters
        ----------
        radian : bool, optional
            If True, returns coordinates in radians; otherwise, degrees.

        Returns
        -------
        tuple of np.ndarray
            Latitude, longitude, and pixel area.
        """
        ny, nx = self.shape
        R = radius.to_value("cm")
        y, x = np.mgrid[0:ny, 0:nx]

        world = self.wcs.pixel_to_world(x, y)
        lat, lon = world.lat, world.lon

        lat, lon = (lat.radian, lon.radian) if radian else (lat.deg, lon.deg)
        area = (2.0 * np.pi * R / nx) ** 2 * np.cos(world.lat.radian)
        return lat, lon, area

    def parameters(self, monopole_corr: bool = False) -> dict:
        """
        Computes key magnetic field parameters.

        Parameters
        ----------
        monopole_corr : bool, optional
            If True, applies monopole correction.

        Returns
        -------
        dict
            Dictionary containing:
            - Total absolute flux (ABSFLUX)
            - Polar region strengths (POLARNxx, POLARSxx)
            - Axial dipole moment (ADM)
            - Equatorial dipole moment (EDM)
            - Azimuthal average (AZAVG)
        """
        lat, lon, area = self.pixel_coord(radian=True)
        coslat, sinlat = np.cos(lat), np.sin(lat)
        coslon, sinlon = np.cos(lon), np.sin(lon)

        # Apply monopole correction if enabled
        sftmap = self.data - (self.data * coslat).sum() / coslat.sum() if monopole_corr else self.data

        # Recompute in degrees
        lat, lon, area = self.pixel_coord(radian=False)

        # Compute absolute flux
        flux = sftmap * area
        total_flux = np.sum(np.abs(flux))

        # Compute polar field strengths
        polar_strengths = {
            f"POLARN{latlim}": flux[lat > latlim].sum() / area[lat > latlim].sum()
            for latlim in [55, 60, 65, 70]
        }
        polar_strengths.update({
            f"POLARS{latlim}": flux[lat < -latlim].sum() / area[lat < -latlim].sum()
            for latlim in [55, 60, 65, 70]
        })

        # Compute azimuthal average
        azavg = self.data.mean(axis=1)

        # Compute dipole moments
        ny, nx = self.shape
        dlat, dlon = np.pi / ny, 2.0 * np.pi / nx

        axial_dipole = (3.0 / (4.0 * np.pi)) * np.sum(sftmap * coslat * sinlat) * dlat * dlon
        equatorial_dipole = np.hypot(
            (3.0 / (4.0 * np.pi)) * np.sum(sftmap * sinlat ** 2 * coslon) * dlat * dlon,
            (3.0 / (4.0 * np.pi)) * np.sum(sftmap * sinlat ** 2 * sinlon) * dlat * dlon,
        )

        return {
            "FILENAME": str(self.filename),
            **dict(self.header),
            **polar_strengths,
            "EDM": equatorial_dipole,
            "ADM": axial_dipole,
            "AZAVG": azavg,
            "ABSFLUX": total_flux,
        }

    def __repr__(self) -> str:
        """Returns a string representation of the object."""
        return f"SFTMap(Model=AFTv2, Data={self.shape}, Time={self.timestamp})"


    def plot(self, draw_mask: bool = False, draw_limb: bool = True, draw_grid: bool = False,
             central_meridian: bool = False, save: bool = False, outpath: str = None,
             verbose: bool = True, **kwargs) -> None:
        """
        Plots the AFT map with optional overlays including limb, grid, mask, and central meridian.

        Parameters
        ----------
        draw_mask : bool, optional
            If True, draws the mask overlay on the plot (default is False).
        draw_limb : bool, optional
            If True, draws the solar limb (default is True).
        draw_grid : bool, optional
            If True, overlays a coordinate grid (default is False).
        central_meridian : bool, optional
            If True, draws the central meridian if available in header (default is False).
        save : bool, optional
            If True, saves the plot to a file instead of displaying it (default is False).
        outpath : str, optional
            Path to save the plot; if None, uses default filename structure.
        verbose : bool, optional
            If True, prints the save location when `save=True` (default is True).
        **kwargs : dict
            Additional keyword arguments for SunPy’s `plot` function.

        Returns
        -------
        None
            Displays or saves the generated AFT map plot.
        """
        smap = self.sunpymap

        # Initialize the figure
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(projection=smap)

        # Plot the SunPy map
        smap.plot(axes=ax, **kwargs)

        # Configure longitude tick labels
        lon = ax.coords["crln"]
        lon.set_ticks(np.arange(30, 359, 60) * u.deg)
        lon.set_major_formatter("d")

        # Optional overlays
        if draw_limb:
            smap.draw_limb(axes=ax)
        if draw_grid:
            smap.draw_grid(axes=ax)
        if draw_mask and self.mask is not None:
            ax.contour(self.mask, levels=[0], colors="orange", linewidths=1, linestyles="dashed")
        if central_meridian and "CRLN_OBS" in self.header:
            num_points = 100
            constant_lon = SkyCoord(
                self.header["CRLN_OBS"] * u.deg,
                np.linspace(-85, 85, num_points) * u.deg,
                frame=smap.coordinate_frame
            )
            ax.plot_coord(constant_lon, color="white", linewidth=2, linestyle=":")

        # Handle saving or displaying the plot
        if save:
            if outpath is None:
                filepath = smap.date.strftime("AFTmap_%Y%m%d_%H%M%S.png")
            else:
                filepath = Path(outpath).resolve() / smap.date.strftime("%Y/%m/AFTmap_%Y%m%d_%H%M%S.png")
                filepath.parent.mkdir(parents=True, exist_ok=True)

            plt.savefig(filepath, dpi=120, bbox_inches="tight")
            plt.close()

            if verbose:
                print(f"SFTMap saved to {filepath}")
        else:
            plt.show()
