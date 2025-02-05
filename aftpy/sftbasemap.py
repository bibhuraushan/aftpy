"""
This the AFTmap Python Class which reads the AFT map data and
provides various parameters for it. It also help to visualize
AFTmap.
"""
__all__ = ['SFTMap']

from typing import Tuple, Any
from astropy.io.fits import Header
from docutils.io import InputError
from sunpy.util.metadata import MetaDict
import h5py as hdf
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from astropy.time import Time
import astropy.units as u
from astropy.io import fits
import sunpy.sun.constants as constants
from pathlib import Path
from sunpy.map.header_helper import make_heliographic_header
from sunpy.coordinates import get_earth
from numpy import ndarray, dtype
from urllib.parse import urlparse
from io import BytesIO
import requests
import os
from .utilities import detect_sftmodel
from .io.io import read_aftv2

this_directory = Path(__file__).parent


class SFTMap:
    """
    A python class to work with individual AFTmap file.

    Attributes
    ----------
    nx:int
        Number of pixels in the x-direction.
    ny: int
        Number of pixels in the y-direction.
    dlat:float
        Latitude step size.
    dlon: float
        Latitude and longitude step size.
    latd: ndarray
        Latitude grid in degrees.
    lond: ndarray
        Longitude grid in degrees.
    latr: ndarray
        Latitude grid in radian.
    lonr: ndarray
        Longitude grid in radian.
    map_list: list, optional
        List to stored information in AFTmap file.
    contents_info : dict, optional
        Dictionary containing information about the AFTmap file.
    fileext: dict, optional
        Dictionary containing extension for various AFTmap file format.

    Parameters
    ----------
    filepath: str
        Name or url of the file to be loaded or read.
    filetype: str, optional
        Type of the file to load. Default is "h5"
    date_fmt:str, optional
        Date format of the file. Default is "AFTmap_%Y%m%d_%H%M.h5"
    timestamp: dt.datetime, optional
        Timestamp of the file if filetype is "hipft". Default is none otherwise.

    Methods
    ---------
    polarfield(self, monopole_corr: bool = True, latlim: float = 60, **kwargs) -> tuple:
        Calculate the polar field of an AFTmap.

    dipole(self, monopole_corr: object = True) -> tuple:
        Calculate the dipole of an AFTmap.

    convert(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True, **kwargs) -> None:
        A conversion tool for AFTmap.
    plot(self, show_mask: bool = True, save: bool = False, outpath: str = None) -> tuple:
        Plot AFTmap data with color bar and timestamp.
    """

    # nx = 1024
    # ny = 512
    # dlat = np.pi / ny
    # dlon = 2.0 * np.pi / nx
    # latd = np.tile(np.arange(0.5, ny) * 180.0 / ny - 90.0, (nx, 1)).T
    # lond = np.tile(np.arange(0.5, nx) * 360.0 / nx, (ny, 1))
    # latr = np.deg2rad(latd)
    # lonr = np.deg2rad(lond)
    # contents_info = {"aftmap": "AFT Baseline Map.",
    #                  "mask": "Region of Data Assimilation.",
    #                  "vlat": "Theta Component of flows at the surface.",
    #                  "vlon": "Phi Component of flows at the surface.",
    #                  "magmap": "Assimilated magnetogram in Carrington Grid."}
    # fileext = {"aftmap": "h5", "oldaft": "dat", "hipft": "h5"}

    def __init__(self, data:[np.ndarray|str], header=None,
                 mask=None, magmap=None, **kwargs) -> None:
        """Initialize an AFTmap object.
        """
        if isinstance(data,str) and header is None:
            _sftmodel = detect_sftmodel(data)
            _sftdata = read_aftv2(data)
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
    def data(self):
        """

        Returns
        -------

        """
        return self._data
    @data.setter
    def data(self, data):
        if (data is None) and (not hasattr(data, "shape")):
            raise InputError("Data must be provided.")
        self._data = data

    @property
    def header(self):
        """

        Returns
        -------

        """
        return self._header
    @header.setter
    def header(self, header):
        """

        Parameters
        ----------
        header
        """
        if (header is None) or not isinstance(header, (MetaDict, Header, dict)):
            raise InputError("Header or MetaDict are required.")
        self._header = header

    @property
    def mask(self):
        """

        Returns
        -------

        """
        return self._mask
    @mask.setter
    def mask(self, mask):
        if (mask is not None) and (not hasattr(mask, "shape")):
            raise InputError("Mask must be 2D array.")
        self._mask = mask

    def __getattr__(self, item):
        return self.header.get(item.lower(), None)

    @property
    def timestamp(self):
        """

        Returns
        -------

        """
        timestamp = (self.header.get("date_obs", None) or
                     self.header.get("date_obs", None) or
                     self.header.get("date_obs", None) or None)
        return timestamp




    def __repr__(self):
        """String represenation of the AFTMap object.
        """
        return f"SFTMap(Model=AFTv2, Data={self.data.shape} Time={self.timestamp})"


    # ============================================================
    # Field parameters from AFT maps;
    # Polar field and dipole moments
    # ============================================================

    # def polarfield(self, monopole_corr: bool = True, latlim: float = 60, **kwargs) -> tuple:
    #     """
    #     A helper function which returns the the polar fields for nortjern hemisphere and southern
    #     hemisphere in the given latitude limits.
    #
    #     Parameters
    #     ----------
    #     monopole_corr:
    #         Whether to apply a monopole correction to the AFT map or not. Defaults to True.
    #     latlim:
    #         Latitude range in which polar field will be calculated. Defaults to 60.
    #     kwargs:
    #         Additional keyword.
    #
    #     Returns
    #     -------
    #     pf:
    #         2 element tuple with polar field data.
    #
    #     """
    #     coslat = np.cos(self.latr)
    #     _aftmap = self.aftmap
    #     if monopole_corr:
    #         _mp = (_aftmap * coslat).sum() / coslat.sum()
    #         _aftmap = _aftmap - _mp
    #
    #     _aftmap = _aftmap.mean(axis=1)
    #     latstrip = self.latd[:, 0]
    #     indN = np.where(latstrip >= latlim)
    #     indS = np.where(latstrip <= -latlim)
    #     latstrip = np.deg2rad(latstrip)
    #     pfN = (_aftmap[indN] * np.cos(latstrip[indN])).sum() / np.cos(latstrip[indN]).sum()
    #     pfS = (_aftmap[indS] * np.cos(latstrip[indS])).sum() / np.cos(latstrip[indS]).sum()
    #     return pfN, pfS
    #
    # def dipole(self, monopole_corr: object = True) -> tuple:
    #     """
    #     Function to calculate the axial and equatorial dipole moments of
    #     the AFTMap object.
    #     Parameters
    #     ----------
    #     monopole_corr: bool
    #         Whether to correct monopole correction or not. Defaults to True.
    #
    #     Returns
    #     -------
    #     dipole: tuple
    #         A tuple containing the axial and equatorial dipole moments.
    #
    #     """
    #     coslat = np.cos(self.latr)
    #     sinlat = np.sin(self.latr)
    #     coslon = np.cos(self.lonr)
    #     sinlon = np.sin(self.lonr)
    #     _aftmap = self.aftmap
    #     if monopole_corr:
    #         _mp = (_aftmap * coslat).sum() / coslat.sum()
    #         _aftmap = _aftmap - _mp
    #
    #     # Axial dipole
    #     _adipole = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * coslat * sinlat) * self.dlat * self.dlon
    #
    #     # Equtorial Dipole
    #     _edipolex = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * sinlat * sinlat * coslon) * self.dlat * self.dlon
    #     _edipoley = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * sinlat * sinlat * sinlon) * self.dlat * self.dlon
    #     _edipole = np.hypot(_edipolex, _edipoley)
    #     return _adipole, _edipole
    #
    # @property
    # def area(self):
    #     """
    #     Calculate the area of each pixel grids.
    #     Returns
    #     -------
    #     area : ndarray
    #         The area of each pixel grids.
    #     """
    #     R = constants.radius.value * 100.0
    #     area = (2.0 * np.pi * R / 1024.0) * (2.0 * np.pi * R / 1024.0) * np.cos(self.latr)
    #     return area
    #
    # @property
    # def flux(self):
    #     """
    #     Flux of each pixel grids.
    #     Returns
    #     -------
    #     flux : ndarray
    #         The flux of each pixel grids.
    #     """
    #     return self.aftmap * self.area
    #
    # @property
    # def pfN(self):
    #     """
    #     Polar field for Northern Hemisphere.
    #     Returns
    #     -------
    #     pfN: float
    #         The polar field for Northern Hemisphere.
    #     """
    #     return self.polarfield()[0]
    #
    # @property
    # def pfS(self):
    #     """
    #     Polar field for Southern Hemisphere for AFT map data.
    #     Returns
    #     -------
    #     pfS: float
    #         The polar field for Southern Hemisphere.
    #     """
    #     return self.polarfield()[1]
    #
    # @property
    # def ADP(self):
    #     """
    #     Axial Dipole Moment of the AFTMap data.
    #     Returns
    #     -------
    #     adp: float
    #         ADP of the Sun.
    #     """
    #     return self.dipole()[0]
    #
    # @property
    # def EDP(self):
    #     """
    #     Equitorial Dipole Moment (EDP) of the AFT map data.
    #     Returns
    #     -------
    #     EDP: float
    #         EDP of the Sun.
    #     """
    #     return self.dipole()[1]
    #
    # # ============================================================
    # # Conversion of AFT map in various formats
    # # ============================================================
    #
    # def convert(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True, **kwargs) -> None:
    #     """
    #     Function to convert all the files in the given directory and save them in
    #     outpath directory with the proper directory structure.
    #
    #     Parameters
    #     ----------
    #     convert_to: str, optionals
    #         The type of converted files.
    #     outpath: str, optionals
    #         The path to save the converted files.
    #     verbose:bool, optionals
    #
    #     Returns
    #     -------
    #     None
    #     """
    #     if self.filetype == "aftmap":
    #         header = self.metadata
    #         data = self.aftmap
    #         if convert_to == "fits":
    #             hdu = fits.PrimaryHDU(data, self.header)
    #             hdu.writeto(outpath, overwrite=True)
    #             if verbose:
    #                 print(f"Output written to {outpath}.")
    #         elif convert_to == "png":
    #             show_mask = kwargs["show_mask"] if "show_mask" in kwargs else False
    #             fig, ax = self.plot(show_mask=show_mask, save=True, outpath=outpath)
    #         else:
    #             print("No implemented Yet.")
    #
    # # ============================================================
    # # Visulaisation of AFT maps;
    # # ============================================================
    # def plot(self, show_mask: bool = True, save: bool = False, outpath: str = None) -> tuple:
    #     """
    #     Display or save the AFT map visualization.
    #
    #     Parameters
    #     ----------
    #     show_mask: bool, optional
    #         Whether to show binary amsk or not.
    #     save :bool, optional
    #         Whether to save the visualization as an image. Defaults to False.
    #     outpath: str
    #         The path to save the.
    #
    #     Returns
    #     ----------
    #     fig: matplotlib.figure.Figure
    #         The Figure object.
    #     ax: matplotlib.axes._axes.Axes
    #         The AxesSubplot object.
    #     """
    #     fig, ax = plt.subplots(figsize=(7, 3.5))
    #     ax.grid(linestyle="--", color="white", alpha=0.3)
    #     plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    #     ax.set_xlabel("Longitude [degree]")
    #     ax.set_ylabel("Latitude [degree]")
    #     _extent = (0, 360, -90, 90)
    #
    #     bmap = self.aftmap
    #
    #     icm = ax.imshow(bmap, origin="lower", extent=_extent,
    #                     cmap="hmimag", vmax=1000, vmin=-1000)
    #
    #     axc = fig.add_axes([0.92, 0.1, 0.025, 0.8])
    #     cb = plt.colorbar(icm, cax=axc, label=r"$B$ [G]")
    #     axc.grid(False)
    #     axc.minorticks_off()
    #     _ticks = [-750, -500, -250, 0, 250, 500, 750]
    #     axc.set_yticks(_ticks, _ticks, rotation="vertical", va="center")
    #     axc.tick_params(direction="out")
    #     if show_mask & (self.filetype == "aftmap"):
    #         ax.contour(self.mask, extent=_extent,
    #                    colors="white", linewidths=1,
    #                    linestyles=":")
    #     ax.set_title(self.time)
    #
    #     if save:
    #         if outpath is not None:
    #             plt.savefig(f"{outpath}/" + "".join(self.name.split(".")[:-1]) + ".png")
    #         else:
    #             plt.savefig("".join(self.name.split(".")[:-1]) + ".png")
    #         plt.close()
    #     else:
    #         plt.show()
    #     return fig, ax
    #
    # # ============================================================
    # # Special functions for AFT maps;
    # # ============================================================
    # @header.setter
    # def header(self, value):
    #     self._header = value
    #
    # @magmap.setter
    # def magmap(self, value):
    #     self._magmap = value
