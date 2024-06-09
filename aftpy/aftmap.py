from typing import Tuple, Any
import h5py as hdf
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import pandas as pd
from astropy.time import Time
import astropy.units as u
from astropy.io import fits
import time
import os
import sunpy.sun.constants as constants
import tqdm as tq
import sunpy.visualization.colormaps
from pathlib import Path
from sunpy.map.header_helper import make_heliographic_header
from sunpy.coordinates import get_earth
from sunpy.coordinates.sun import carrington_rotation_number
from .utilities import file_search, h52df, df2h5
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import cpu_count
from .visulalization import Visulalization
from numpy import ndarray, dtype

this_directory = Path(__file__).parent

__all__ = ['AFTmap', 'AFTmaps']


class AFTmap:
    """
    A class to work with AFT maps.

    Attributes:
    - nx (int): Number of pixels in the x-direction.
    - ny (int): Number of pixels in the y-direction.
    - dlat (float): Latitude step size.
    - dlatlon (float): Latitude and longitude step size.

    Methods:
    - aftmap (property): Get the AFT map.
    - mask (property): Get the mask.
    - vv (property): Get the vlat and vlon.
    - metadata (property): Get the metadata.
    - info (property): Get additional information.
    - __str__(): Return a string representation of the object.
    - view(): Display or save the AFT map visualization.
    """
    nx = 1024
    ny = 512
    dlat = np.pi / ny
    dlon = 2.0 * np.pi / nx
    latd = np.tile(np.arange(0.5, ny) * 180.0 / ny - 90.0, (nx, 1)).T
    lond = np.tile(np.arange(0.5, nx) * 360.0 / nx, (ny, 1))
    latr = np.deg2rad(latd)
    lonr = np.deg2rad(lond)
    contents_info = {"aftmap": "AFT Baseline Map.",
                     "mask": "Region of Data Assimilation.",
                     "vlat": "Theta Component of flows at the surface.",
                     "vlon": "Phi Component of flows at the surface.",
                     "magmap": "Assimilated magnetogram in Carrington Grid."}
    fileext = {"aftmap": "h5", "oldaft": "dat", "hipft": "h5"}

    def __init__(self, file: str, filetype: str = "aftmap",
                 date_fmt: str = "AFTmap_%Y%m%d_%H%M.h5",
                 timestamp: dt.datetime = None):
        """
        Initialize an AFTmap object.

        Parameters
        ----------
        file: str
        Name of the to be loaded or read.
        filetype: str, optional
        Type of the file to load. Default is "h5"
        date_fmt:str, optional
        Date format of the file. Default is "AFTmap_%Y%m%d_%H%M.h5"
        timestamp: dt.datetime, optional
        Timestamp of the file if filetype is "hipft". Default is dt.datetime
        """
        self.file = file
        self.name = file.split("/")[-1]
        self.filetype = filetype
        self.date_fmt = date_fmt
        self.map_list = None
        self.timestamp = timestamp

        if self.filetype == "aftmap":
            with hdf.File(self.file) as fl:
                self.map_list = [_key for _key in fl["maps"].keys()]
        else:
            self.map_list = ["aftmap"]

        if (filetype == "hipft") & (self.timestamp is None):
            print("Timestamp not provided. Assuming today as the timestamp.")
            self.timestamp = Time(dt.datetime.today()).fits
        elif timestamp is not None:
            self.timestamp = Time(timestamp).fits

    @property
    def contents(self) -> list:
        """
        Returns the contents of the AFTmap file as list.

        Returns
        -------
        maplist: list
        Contents of the AFTmap file.

        """
        return self.map_list

    # ============================================================
    # Data attributes associated with file
    # ============================================================
    @property
    def time(self) -> str:
        if self.filetype == "aftmap":
            return self.metadata["map_date"]
        elif self.filetype == "oldaft":
            _time = Time(dt.datetime.strptime(
                self.name, self.date_fmt), scale="utc")
        elif self.filetype == "hipft":
            _time = self.timestamp
            return _time

    @property
    def ymd(self) -> tuple:
        year, month, day = self.time.split("-")[0:3]
        return year, month, day[0:2]

    @property
    def aftmap(self) -> np.ndarray:
        """
        Get the AFT map data.

        Returns
        -------
        aftmap: ndarray
        The AFT map file.
        """
        if self.filetype == "aftmap":
            # For HDF file
            with hdf.File(self.file) as fl:
                bmap = np.array(fl["maps/aftmap"])
        elif self.filetype == "oldaft":
            # For old dat file
            _data = np.fromfile(self.file, dtype=np.float32)
            if _data.size == 512 * 1024:
                bmap = _data.reshape(512, -1)
            else:
                bmap = _data[0:int(512 * 1024)].reshape(512, -1)
        elif self.filetype == "hipft":
            with hdf.File(self.file) as fl:
                bmap = np.array(fl["Data"])
        else:
            bmap = None
        return bmap

    @property
    def mask(self) -> np.ndarray:
        """
        Get the mask data.

        Returns
        -------
        mask: ndarray
        The region where data assimilation has been performed.
        """
        if "mask" in self.map_list:
            with hdf.File(self.file) as fl:
                mask = np.array(fl["maps/mask"])
        elif self.filetype == "dat":
            raise Warning("No Mask File found in dat file.")
        else:
            mask = None
        return mask

    @property
    def magmap(self) -> np.ndarray:
        """
        Get the hmi assimilated data if available.

        Returns
        -------
        magmap: ndarray
        The hmi assimilated data if available.
        """
        if "magmap" in self.map_list:
            with hdf.File(self.file) as fl:
                mask = np.array(fl["maps/magmap"])
        elif self.filetype == "dat":
            raise Warning("No Mask File found in dat file.")
        else:
            mask = None
        return mask

    @property
    def vlat(self) -> ndarray[Any, dtype[Any]] | None:
        """
        Get vlat data.

        Returns:
        - vlat (numpy.ndarray): Vlat data.
        """
        if "vlat" in self.map_list:
            with hdf.File(self.file) as fl:
                vlat = np.array(fl["maps/vlat"])
        else:
            vlat = None
        return vlat

    @property
    def vlon(self) -> ndarray[Any, dtype[Any]] | None:
        """
        Get vlon data.

        Returns
        -------
        vlon: ndarray[Any, dtype[Any]]
        Vlon data.
        """
        if "vlon" in self.map_list:
            with hdf.File(self.file) as fl:
                vlon = np.array(fl["maps/vlon"])
        else:
            vlon = None
        return vlon

    # ============================================================
    # Metadata and header from AFT maps
    # ============================================================
    @property
    def metadata(self) -> dict[str, Any]:
        """
        Get metadata.

        Returns:
        - header (dict): Metadata.
        """
        header = {}
        if self.filetype == "aftmap":
            with hdf.File(self.file) as fl:
                for _key in fl["header"].keys():
                    _val = list(fl["header"][_key])[0]
                    if isinstance(_val, bytes):
                        _val = _val.decode("utf-8")
                    header[_key] = _val
        else:
            header = None
        return header

    @property
    def header(self) -> fits.header.Header:
        _headeraft = self.metadata
        date = self.time
        observer = get_earth(date)
        _header = make_heliographic_header(date, observer, (512, 1024),
                                           frame='carrington',
                                           projection_code="CAR",
                                           map_center_longitude=180 * u.deg)
        _header['hgln_obs'] = self.metadata["crln_obs"]
        _header['hglt_obs'] = self.metadata["crlt_obs"]
        header = fits.header.Header(_header)
        header.update(_headeraft)
        return header

    @property
    def info(self) -> dict[str, Any]:
        """
        Get additional information.

        Returns:
        - info (dict): Additional information.
        """
        info = {}
        if self.filetype == "aftmap":
            with hdf.File(self.file) as fl:
                for _key in fl["header"].keys():
                    try:
                        _info = list(fl["header"][_key].attrs["descr"])[0].decode("utf-8")
                    except:
                        _info = "None"
                    info[_key] = _info
        else:
            info = None
        return info

    # ============================================================
    # Field parameters from AFT maps;
    # Polar field and dipole moments
    # ============================================================

    def polarfield(self, monopole_corr: bool = True, latlim: float = 60, **kwargs):
        """
        A helper function which returns the the polar fields for nortjern hemisphere and southern
        hemisphere in the given latitude limits.

        Parameters
        ----------
        monopole_corr: bool
        Whether to apply a monopole correction to the AFT map or not. Defaults to True.
        latlim: float
        Latitude range in which polar field will be calculated. Defaults to 60.
        kwargs: dict
        Additional keyword.

        Returns
        -------
        pf: list
        2 dimensional array of type list.

        """
        coslat = np.cos(self.latr)
        _aftmap = self.aftmap
        if monopole_corr:
            _mp = (_aftmap * coslat).sum() / coslat.sum()
            _aftmap = _aftmap - _mp

        _aftmap = _aftmap.mean(axis=1)
        latstrip = self.latd[:, 0]
        indN = np.where(latstrip >= latlim)
        indS = np.where(latstrip <= -latlim)
        latstrip = np.deg2rad(latstrip)
        pfN = (_aftmap[indN] * np.cos(latstrip[indN])).sum() / np.cos(latstrip[indN]).sum()
        pfS = (_aftmap[indS] * np.cos(latstrip[indS])).sum() / np.cos(latstrip[indS]).sum()
        return [pfN, pfS]

    def dipole(self, monopole_corr: object = True) -> tuple:
        """
        Function to calculate the axial and equatorial dipole moments of
        the AFTMap object.
        Parameters
        ----------
        monopole_corr: bool
        Whether to correct monopole correction or not. Defaults to True.

        Returns
        -------
        dipole: tuple
        A tuple containing the axial and equatorial dipole moments.

        """
        coslat = np.cos(self.latr)
        sinlat = np.sin(self.latr)
        coslon = np.cos(self.lonr)
        sinlon = np.sin(self.lonr)
        _aftmap = self.aftmap
        if monopole_corr:
            _mp = (_aftmap * coslat).sum() / coslat.sum()
            _aftmap = _aftmap - _mp

        # Axial dipole
        _adipole = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * coslat * sinlat) * self.dlat * self.dlon

        # Equtorial Dipole
        _edipolex = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * sinlat * sinlat * coslon) * self.dlat * self.dlon
        _edipoley = (3.0 / (4.0 * np.pi)) * np.sum(_aftmap * sinlat * sinlat * sinlon) * self.dlat * self.dlon
        _edipole = np.hypot(_edipolex, _edipoley)
        return _adipole, _edipole

    @property
    def area(self):
        """
        Calculate the area of each pixel grids.
        Returns
        -------
        area : ndarray
        The area of each pixel grids.
        """
        R = constants.radius.value * 100.0
        area = (2.0 * np.pi * R / 1024.0) * (2.0 * np.pi * R / 1024.0) * np.cos(self.latr)
        return area

    @property
    def flux(self):
        """
        Flux of each pixel grids.
        Returns
        -------
        flux : ndarray
        The flux of each pixel grids.
        """
        return self.aftmap * self.area

    @property
    def pfN(self):
        """
        Polar field for Northern Hemisphere.
        Returns
        -------
        pfN: float
        The polar field for Northern Hemisphere.
        """
        return self.polarfield()[0]

    @property
    def pfS(self):
        """
        Polar field for Southern Hemisphere.
        Returns
        -------
        pfS: float
        The polar field for Southern Hemisphere.
        """
        return self.polarfield()[1]

    @property
    def ADP(self):
        """
        ADP of the Sun.
        Returns
        -------
        adp: float
        ADP of the Sun.
        """
        return self.dipole()[0]

    @property
    def EDP(self):
        """
        EDP of the Sun.
        Returns
        -------
        EDP: float
        EDP of the Sun.
        """
        return self.dipole()[1]

    # ============================================================
    # Conversion of AFT map in various formats
    # ============================================================

    def convert(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True, **kwargs) -> None:
        """
        Function to convert all the files in the given directory and save them in
        outpath directory with the proper directory structure.

        Parameters
        ----------
        convert_to: str, optionals
        The type of converted files.
        outpath: str, optionals
        The path to save the converted files.
        verbose:bool, optionals

        Returns
        -------
        object
        Whether to show progress.
        """
        if self.filetype == "aftmap":
            header = self.metadata
            data = self.aftmap
            if convert_to == "fits":
                hdu = fits.PrimaryHDU(data, self.header)
                hdu.writeto(outpath, overwrite=True)
                if verbose:
                    print(f"Output written to {outpath}.")
            elif convert_to == "png":
                show_mask = kwargs["show_mask"] if "show_mask" in kwargs else False
                fig, ax = self.plot(show_mask=show_mask, save=True, outpath=outpath)
            else:
                print("No implemented Yet.")

    # ============================================================
    # Visulaisation of AFT maps;
    # ============================================================
    def plot(self, show_mask: bool = True, save: bool = False, outpath=None) -> tuple:
        """
        Display or save the AFT map visualization.

        Args:
        - para (str, optional): Parameter to display ("aftmap" or "mask"). Defaults to "aftmap".
        - save (bool, optional): Whether to save the visualization as an image. Defaults to False.

        Returns:
        - fig (matplotlib.figure.Figure): The Figure object.
        - ax (matplotlib.axes._axes.Axes): The AxesSubplot object.

        Parameters
        ----------
        outpath: str
        The path to save the.
        show_mask: bool, optional
        Whether to show the mask of the AFT map. Defaults to True
        save: bool, optional
        Whether to save the visualization as an image. Defaults to False.

        Returns
        -------
        fig, object:
        The Figure object.
        ax: matplotlib.axes._axes.AxesSubplot
        The AxesSubplot object.
        """
        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.grid(linestyle="--", color="white", alpha=0.3)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        ax.set_xlabel("Longitude [degree]")
        ax.set_ylabel("Latitude [degree]")
        _extent = (0, 360, -90, 90)

        bmap = self.aftmap

        icm = ax.imshow(bmap, origin="lower", extent=_extent,
                        cmap="hmimag", vmax=1000, vmin=-1000)

        axc = fig.add_axes([0.92, 0.1, 0.025, 0.8])
        cb = plt.colorbar(icm, cax=axc, label=r"$B$ [G]")
        axc.grid(False)
        axc.minorticks_off()
        _ticks = [-750, -500, -250, 0, 250, 500, 750]
        axc.set_yticks(_ticks, _ticks, rotation="vertical", va="center")
        axc.tick_params(direction="out")
        if show_mask & (self.filetype == "aftmap"):
            ax.contour(self.mask, extent=_extent,
                       colors="white", linewidths=1,
                       linestyles=":")
        ax.set_title(self.time)

        if save:
            if outpath is not None:
                plt.savefig(f"{outpath}/" + "".join(self.name.split(".")[:-1]) + ".png")
            else:
                plt.savefig("".join(self.name.split(".")[:-1]) + ".png")
            plt.close()
        else:
            plt.show()
        return fig, ax

    # ============================================================
    # Special functions for AFT maps;
    # ============================================================

    def __repr__(self):
        """
        String represenation of the AFTMap object.
        Returns
        -------
        str: str
        String representation of the AFTMap object.
        """
        return f"AFTMap(time={self.time}, filtype={self.filetype})"

    def __str__(self) -> str:
        """
        Get a string representation of the object.

        Returns
        -------
        str
        String representation of the object.
        """
        dct = self.metadata
        dct1 = self.info
        print(" " * 30 + f"{self.name}" + " " * 30)
        if self.filetype == "aftmap":
            print("-" * 85)
            print(" " * 36 + f"MAP CONTENTS" + " " * 37)
            print("-" * 85)
            for _key in self.map_list:
                print("%-10s: %-30s - %-30s" % (_key.upper(), "Array", self.contents_info[_key]))
            print("-" * 85)
            print(" " * 39 + f"HEADER" + " " * 40)
            print("-" * 85)
            for _key in self.info.keys():
                print("%-10s: %-30s - %-30s" % (_key.upper(), dct[_key], dct1[_key]))
        print("-" * 85)
        print(" " * 39 + f"AFT MAP" + " " * 40)
        print("-" * 85)
        self.plot()
        return "-" * 85


class AFTmaps:
    """
        A class for loading AFT maps from a directory.

        Attributes:
        - path (str): The path to the directory containing AFT map files.
        - filetype (str): The file extension of the AFT map files (e.g., "h5").
        - date_fmt (str): The date format string used to parse timestamps from filenames.
        - filelist (list): A list of file paths to AFT map files in the specified directory.
        - filenames (numpy.ndarray): An array of filenames extracted from the filelist.

        Methods:
        - __init__(path=".", filetype="h5", date_fmt="AFTmap_%Y%m%d_%H%M.h5", verbose=True): Initializes the AFTload object.
        - stats(): Prints statistics about the loaded AFT map files.
    """

    def __init__(self, path: str | list | tuple = ".", filetype: str = "aftmap",
                 date_fmt: str = "AFTmap_%Y%m%d_%H%M.h5",
                 verbose: bool = True, hipft_prop: dict = None):
        """
        Initilaizing class function.

        Parameters
        ----------
        path: str, optional
        The path to the AFT map file. Defaults to the current working directory.
        filetype: str, optional
        The file extension of the AFT map. Default to the "h5"
        date_fmt:str, optional
        The date format string. Defaults to the "AFTmap_%Y%m%d_%H%M.h5"
        verbose:bool, optional
        Whether to print or not. Defaults to the True
        hipft_prop: dict, optional
        A dictionary for hipft properties, if filetype is "hipft", the hipft
        """
        self.path = path
        self.filetype = filetype
        self.date_fmt = date_fmt
        self._fileext = {"aftmap": "h5", "oldaft": "dat", "hipft": "h5"}

        # Whether it is a list of paths or one path
        if isinstance(path, (list, tuple)):
            filelist = np.array([])
            for _path in path:
                _filelist = file_search(_path, fileext=self._fileext[filetype])
                filelist = np.append(filelist, _filelist)
            self.filelist = filelist
        else:
            self.filelist = file_search(path, fileext=self._fileext[filetype])

        self.filenames = np.array([os.path.basename(f) for f in self.filelist])
        self.hipft_prop = hipft_prop
        self.counts = len(self.filenames)
        if (filetype == "hipft") & (hipft_prop is None):
            raise Exception("HIPFT Time information not provided")

        # Carrington Parameters
        self._crn = np.ceil(carrington_rotation_number(self.timestamps)).astype(int)
        self._CRN = np.unique(self._crn).astype(int)

        if verbose:
            self.stats()

    def __len__(self) -> int:
        """
        Length of the object.
        Returns
        -------
        n: int
        Length of the object.
        """
        return len(self.filelist)

    def __repr__(self):
        return f"<AFTMaps(Path: {self.path}, Type: {self.filetype}, Total: {self.counts})>"

    def __str__(self):
        self.stats()

    @property
    def crn(self):
        return self._CRN

    @property
    def aftmaps(self):
        """
        An itterator containing all the AFTmaps associated.
        """
        for _file in self.filelist:
            yield AFTmap(_file, date_fmt=self.date_fmt, filetype=self.filetype)

    @property
    def timestamps(self):
        """
        To get timestamps for all files in the AFT map directory.

        Returns
        -------
        object
        List of timestamps corresponding to each AFT map in the directory.

        """
        _timestamp = []
        if self.filetype != "hipft":
            for _file in self.filelist:
                _timestamp.append(dt.datetime.strptime(os.path.basename(_file), self.date_fmt))
        else:
            _t0 = self.hipft_prop["T0"]
            _dt = self.hipft_prop["dt"]
            _timestamp = [_t0 + dt.timedelta(i * _dt) for i in range(len(self.filelist))]
        return Time(_timestamp)

    def get_filelist(self) -> np.ndarray:
        """
        To get the full paths of the loaded files in the AFTmaps folder.

        Returns
        -------
        np.ndarray
        List of all files in the AFTmaps directory.

        """
        return self.filelist

    def get_filenames(self) -> np.ndarray:
        """
        Returns the names of the files in the AFTMap directory.

        Returns
        -------
        np.ndarray
        list of filenames.

        """
        return np.array(self.filenames)

    def stats(self):
        """
        Show the statistics of the loaded files.
        """
        _stats = {"RootDir": self.path, "FileType": self.filetype, "# Files": len(self.filelist),
                  "T-Initial": self.timestamps[0], "T-End": self.timestamps[-1]}
        for _s in _stats.keys():
            print("%-10s: %-30s" % (_s, _stats.get(_s)))

    def convert_all(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True, **kwargs):
        """
        Convert all loaded AFT map files to the specified format.

        Parameters
        ----------
        convert_to: str, optional
        The output format to AFTMap to fits. Defaults to "fits"
        outpath: str, optional
        The directory in which to save the fits. Defaults to "."
        verbose: bool, optional
        Show progress. Defaults to True
        """
        for _file, i in zip(self.filelist, range(len(self.filelist))):
            mapobj = AFTmap(_file)
            yr, mo, day = mapobj.ymd
            _outpath = os.path.join(outpath, str(yr) + "/" + str(mo))
            _filename = os.path.splitext(os.path.basename(_file))[0] + "." + convert_to
            if not os.path.exists(_outpath):
                os.makedirs(_outpath, exist_ok=True)
            if convert_to == "png":
                mapobj.convert(convert_to=convert_to, verbose=False, outpath=_outpath, **kwargs)
            else:
                _filename = os.path.join(_outpath, _filename)
                mapobj.convert(convert_to=convert_to, verbose=False, outpath=_filename, **kwargs)
            if verbose:
                frac = float(i + 1) / len(self.filelist) * 100
                print(f"({frac:.2f}%) Converting {os.path.basename(_file)} to {os.path.basename(_filename)}", end="\r")

    @staticmethod
    def _generate_para(file):
        """
        A helper function to generate a aft parameters for a given file.
        Parameters
        ----------
        file: str
        The name of the file to generate the aft parameters.

        Returns
        -------
        para: tuple
        A tuple containing the all the AFT parameters for the file.

        """
        mapobj = AFTmap(file)
        tflux = round(np.abs(mapobj.flux).sum(), 3)
        pf = mapobj.polarfield()
        dp = mapobj.dipole()
        para = (mapobj.time, *pf, *dp, tflux)
        return para

    def generate_parameters(self, outfile=None, nthreds=None,
                            verbose=True, monopole_corr=True, use_saved=True):
        """
        Function to generate AFT data for all the files.
        Parameters
        ----------
        use_saved: bool, optional
        Whether to relaod form saved file or not. Defaults to True.
        monopole_corr: bool, optional
        Whether to apply a monopole correction to the data. Default is True
        outfile: str, optional
        The file in which data is to be saved. Defaults to None.
        nthreds: int, optional
        The number of thredds to use for the calculations. Defaults to None.
        verbose: bool, optional
        Show progress. Defaults to True.

        Returns
        -------
        df: pd.DataFrame
        The dataframe containing the AFT parameters for the given files.

        """
        if nthreds is None:
            nthreds = cpu_count() - 1

        save_file = os.path.join(this_directory, ".aftpara.h5")
        save_exist = os.path.isfile(save_file)

        if use_saved & save_exist:
            df = h52df(save_file)
            if (df.Time.iloc[-1] == self.timestamps[-1]) and (self.counts == df.shape[0]):
                return df
            elif verbose:
                print("WARNING: Saved data is outdated, skipping saved data.")

        if outfile is None:
            outfile = dt.datetime.now().strftime("AFTpara_%Y%m%d_%H%M.csv")
        result = []
        tint = time.time()
        bar_format = '{desc}: {percentage:3.2f}%|{bar}{r_bar}'
        with Pool(nthreds) as pool:
            for _para in tq.tqdm(pool.imap(self._generate_para, self.filelist),
                                 bar_format=bar_format, total=self.counts):
                result.append(_para)
        if verbose:
            print(f"Completed in {int(time.time() - tint)} seconds.")
            print(f"AFT parameters are saved to {outfile}.")
        df = pd.DataFrame(result, columns=["Time", "PolarN", "PolarS", "ADM", "EDM", "TotalFlux"])
        df2h5(df, save_file)
        if outfile is not None:
            df.to_csv(outfile, float_format="%.2g", index=False)
        return df

    def get_crmap(self, cr_number):
        if cr_number not in self._CRN:
            raise ValueError("The given Carrington Rotation is Not present in the data.")

        ind = np.where(self._crn == cr_number)
        crmap = 0
        for _fl in self.filelist[ind]:
            aftobj = AFTmap(_fl, filetype=self.filetype, date_fmt=self.date_fmt)
            crmap = crmap + aftobj.aftmap
        crmap = crmap / ind[0].size
        return crmap

    def cravgmap(self, verbose=False, outfile=None):
        save_file = os.path.join(this_directory, ".bflydata.h5")
        save_exist = os.path.isfile(save_file)
        # if use_saved & save_exist:

        # Get carrington Numbers
        crn = self._CRN
        crn0 = crn.min()
        bar_format = '{desc}: {percentage:3.2f}%|{bar}{r_bar}'
        bflymap = np.zeros((512, len(crn)))
        for _crn in tq.tqdm(crn, bar_format=bar_format, total=len(crn)):
            bflymap[:, _crn - crn0] = self.get_crmap(_crn).mean(axis=1)
        return bflymap

    @property
    def visualize(self, **kwargs):
        return Visulalization(self, **kwargs)
