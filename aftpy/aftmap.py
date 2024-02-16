from typing import Tuple, Any

import h5py as hdf
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from astropy.time import Time
from astropy.io import fits
import os
from pathlib import Path

from numpy import ndarray, dtype

this_directory = Path(__file__).parent
import sunpy.visualization.colormaps

try:
    plt.style.use(f'{this_directory}/bkj_style.mplstyle')
except:
    plt.style.use('default')
    print("Using default matplotlib style.")


def file_search(path, fileext: str = None) -> np.ndarray:
    """
       Search for files with the specified file extension in the given directory.

       Args:
       - path (str): The directory path to search for files.
       - fileext (str, optional): The file extension to filter files. Defaults to "" (no filtering).

       Returns:
       - filelist (list): A list of file paths matching the specified criteria.
    """
    files = []
    for dr, _, file in os.walk(path):
        for f in file:
            if fileext is None:
                files.append(os.path.join(dr, f))
            elif f.endswith(fileext):
                files.append(os.path.join(dr, f))
    return np.sort(np.array(files))


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
    dlatlon = 2.0 * np.pi / nx
    contents_info = {"aftmap": "AFT Baseline Map.",
                     "mask": "Region of Data Assimilation.",
                     "vlat": "Theta Component of flows at the surface.",
                     "vlon": "Phi Component of flows at the surface.",
                     "magmap": "Assimilated magnetogram in Carrington Grid."}

    def __init__(self, file, filetype: str = "h5",
                 date_fmt: str = "AFTmap_%Y%m%d_%H%M.h5",
                 timestamp: dt.datetime = None):
        """
        Initialize an AFTmap object.

        Args:
        - file (str): Path to the HDF5 file.
        """
        self.file = file
        self.name = file.split("/")[-1]
        self.filetype = filetype
        self.date_fmt = date_fmt
        self.map_list = None
        self.timestamp = timestamp

        if self.filetype == "h5":
            with hdf.File(self.file) as fl:
                self.map_list = [_key for _key in fl["maps"].keys()]

        if (filetype == "hipft") & (self.timestamp is None):
            print("Timestamp not provided. Assuming today as the timestamp.")
            self.timestamp = Time(dt.datetime.today()).fits
        elif timestamp is not None:
            self.timestamp = Time(timestamp).fits

    @property
    def contents(self) -> list:
        return self.map_list

    @property
    def aftmap(self) -> np.ndarray:
        """
        Get the AFT map data.

        Returns:
        - bmap (numpy.ndarray): AFT map data.
        """
        if self.filetype == "h5":
            # For HDF file
            with hdf.File(self.file) as fl:
                bmap = np.array(fl["maps/aftmap"])
        elif self.filetype == "dat":
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

        Returns:
        - mask (numpy.ndarray): Mask data.
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
    def vv(self) -> tuple[ndarray[Any, dtype[Any]] | None, ndarray[Any, dtype[Any]] | None]:
        """
        Get vlat and vlon data.

        Returns:
        - vlat (numpy.ndarray): Vlat data.
        - vlon (numpy.ndarray): Vlon data.
        """
        if "vlat" in self.map_list:
            with hdf.File(self.file) as fl:
                vlat = np.array(fl["maps/vlat"])
                vlon = np.array(fl["maps/vlon"])
        else:
            vlat, vlon = None, None
        return vlat, vlon

    @property
    def metadata(self) -> dict[str, Any]:
        """
        Get metadata.

        Returns:
        - header (dict): Metadata.
        """
        header = {}
        if self.filetype == "h5":
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
    def time(self) -> str:
        if self.filetype == "h5":
            return self.metadata["map_date"]
        elif self.filetype == "dat":
            _time = Time(dt.datetime.strptime(
                self.name, self.date_fmt))
        elif self.filetype == "hipft":
            _time = self.timestamp
            return _time

    @property
    def ymd(self) -> tuple:
        year, month, day = self.time.split("-")[0:3]
        return year, month, day[0:2]

    @property
    def info(self) -> dict[str, Any]:
        """
        Get additional information.

        Returns:
        - info (dict): Additional information.
        """
        info = {}
        if self.filetype == "h5":
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

    def convert(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True):
        if self.filetype == "h5":
            header = self.metadata
            data = self.aftmap
            if convert_to == "fits":
                hdu = fits.PrimaryHDU(data)
                hdu.header.update(header)
                hdu.writeto(outpath)
                if verbose:
                    print(f"Output written to {outpath}.")
            else:
                print("No implemented Yet.")

    def plot(self, show_mask: bool = True, save: bool = False) -> tuple:
        """
        Display or save the AFT map visualization.

        Args:
        - para (str, optional): Parameter to display ("aftmap" or "mask"). Defaults to "aftmap".
        - save (bool, optional): Whether to save the visualization as an image. Defaults to False.

        Returns:
        - fig (matplotlib.figure.Figure): The Figure object.
        - ax (matplotlib.axes._axes.Axes): The AxesSubplot object.
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
        if show_mask & (self.filetype == "h5"):
            ax.contour(self.mask, extent=_extent,
                       colors="white", linewidths=1,
                       linestyles=":")
        ax.set_title(self.time)

        if save:
            plt.savefig("data/" + "".join(self.name.split(".")[:-1]) + ".png")
            print("data/" + "".join(self.name.split(".")[:-1]) + ".png")
        else:
            plt.show()
        return fig, ax

    def __str__(self) -> str:
        """
        Get a string representation of the object.

        Returns:
        - str: String representation.
        """
        dct = self.metadata
        dct1 = self.info
        print(" " * 30 + f"{self.name}" + " " * 30)
        if self.filetype == "h5":
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


class AFTload:
    """
        A class for loading AFT maps from files.

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

    def __init__(self, path: str = ".", filetype: str = "h5",
                 date_fmt: str = "AFTmap_%Y%m%d_%H%M.h5",
                 verbose: bool = True, hipft_prop: dict = None):
        """
            Initialize the AFTload object.

            Args:
            - path (str, optional): The path to the directory containing AFT map files. Defaults to current directory.
            - filetype (str, optional): The file extension of the AFT map files. Defaults to "h5".
            - date_fmt (str, optional): The date format string used to parse timestamps from filenames.
              Defaults to "AFTmap_%Y%m%d_%H%M.h5".
            - verbose (bool, optional): Whether to print statistics about the loaded AFT map files. Defaults to True.
        """
        self.path = path
        self.filetype = filetype
        self.date_fmt = date_fmt
        if filetype == "hipft":
            _fileext = "h5"
        else:
            _fileext = filetype
        self.filelist = file_search(path, fileext=_fileext)
        self.filenames = np.array([os.path.basename(f) for f in self.filelist])
        self.hipft_prop = hipft_prop
        if (filetype == "hipft") & (hipft_prop is None):
            raise Exception("HIPFT Time information not provided")

        if verbose:
            self.stats()

    @property
    def aftmaps(self):
        for _file in self.filelist:
            yield AFTmap(_file, date_fmt=self.date_fmt, filetype=self.filetype)

    @property
    def timestamps(self):
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
        return self.filelist

    def get_filenames(self) -> np.ndarray:
        return np.array(self.filenames)

    def stats(self):
        _stats = {}
        _stats["RootDir"] = self.path
        _stats["FileType"] = self.filetype
        _stats["# Files"] = len(self.filelist)
        _stats["T-Initial"] = self.timestamps[0]
        _stats["T-End"] = self.timestamps[-1]
        for _s in _stats.keys():
            print("%-10s: %-30s" % (_s, _stats.get(_s)))

    def convert_all(self, convert_to: str = "fits", outpath: str = ".", verbose: bool = True):
        """
        Convert all loaded AFT map files to the specified format.

        Args:
        - convert_to (str, optional): The output format to convert the AFT map files to. Defaults to "fits".
        - outpath (str, optional): The directory path to save the converted files. Defaults to current directory.
        - verbose (bool, optional): Whether to print conversion progress. Defaults to True.
        """
        for _file, i in zip(self.filelist, range(len(self.filelist))):
            mapobj = AFTmap(_file)
            yr, mo, day = mapobj.ymd
            _outpath = os.path.join(outpath, str(yr) + "/" + str(mo))
            if not os.path.exists(_outpath):
                os.makedirs(_outpath, exist_ok=True)
            _filename = os.path.splitext(os.path.basename(_file))[0] + "." + convert_to
            _filename = os.path.join(_outpath, _filename)
            mapobj.convert(convert_to=convert_to, verbose=False, outpath=_filename)
            if verbose:
                frac = float(i + 1) / len(self.filelist) * 100
                print(f"({frac:.2f}%) Converting {os.path.basename(_file)} to {os.path.basename(_filename)}", end="\r")
