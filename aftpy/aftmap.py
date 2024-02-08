import h5py as hdf
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from astropy.time import Time
from astropy.io import fits
import os
import sunpy.visualization.colormaps

try:
    plt.style.use('aftpy/bkj_style.mplstyle')
except:
    plt.style.use('default')
    print("Using default matplotlib style.")


def file_search(path, fileext=None):
    files = []
    for dr, _, file in os.walk(path):
        for f in file:
            if fileext is None:
                files.append(os.path.join(dr, f))
            elif f.endswith(fileext):
                files.append(os.path.join(dr, f))
    return np.sort(files)


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

    def __init__(self, file, filetype="h5", date_fmt="AFTmap_%Y%m%d_%H%M.h5"):
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

        if self.filetype == "h5":
            with hdf.File(self.file) as fl:
                self.map_list = [_key for _key in fl["maps"].keys()]

    @property
    def contents(self):
        return self.map_list

    @property
    def aftmap(self):
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
        else:
            bmap = None
        return bmap

    @property
    def mask(self):
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
    def vv(self):
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
    def metadata(self):
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
    def time(self):
        if self.filetype == "h5":
            return self.metadata["map_date"]
        elif self.filetype == "dat":
            _time = Time(datetime.datetime.strptime(
                self.name, self.date_fmt))
            return _time

    @property
    def info(self):
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

    def convert(self, outfile="fits", outpath=".", verbose=True):
        if self.filetype == "h5":
            header = self.metadata
            data = self.aftmap
            if outfile == "fits":
                hdu = fits.PrimaryHDU(data)
                hdu.header.update(header)
                hdu.writeto(outpath)
                if verbose:
                    print(f"Output written to {outpath}.")
            else:
                print("No implemented Yet.")

    def plot(self, show_mask=True, save=False):
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

    def __str__(self):
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

    def __init__(self, path=".", filetype="h5", date_fmt="AFTmap_%Y%m%d_%H%M.h5", verbose=True):
        self.path = path
        self.filetype = filetype
        self.date_fmt = date_fmt
        self.filelist = file_search(path, fileext=filetype)
        if verbose: self.stats()

    @property
    def aftmaps(self):
        for _file in self.filelist:
            yield AFTmap(_file, date_fmt=self.date_fmt, filetype=self.filetype)

    @property
    def timestamps(self):
        _timestamp = []
        for _file in self.filelist:
            _timestamp.append(dt.datetime.strptime(os.path.basename(_file), self.date_fmt))
        return Time(_timestamp)

    def stats(self):
        _stats = {}
        _stats["RootDir"] = self.path
        _stats["FileType"] = self.filetype
        _stats["# Files"] = len(self.filelist)
        _stats["T-Initial"] = self.timestamps[0]
        _stats["T-End"] = self.timestamps[-1]
        for _s in _stats.keys():
            print("%-10s: %-30s" % (_s, _stats.get(_s)))

    # def convert_all(self, outfile="fits", outpath=".", verbose=True):
    #     for _file in self.filelist:
