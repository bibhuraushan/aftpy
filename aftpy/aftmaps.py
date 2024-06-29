"""
This class will initialize with getting all the list of files
and create aftmap objects. This will alos provide methods to create
butterfly diagram and generate AFT parameters for all the maps in the list.
"""
import h5py as hdf
import numpy as np
import datetime as dt
import pandas as pd
import time
import os
import tqdm as tq
from pathlib import Path
from sunpy.coordinates.sun import carrington_rotation_number
from .utilities import file_search, h52df, df2h5
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import cpu_count
from .visulalization import Visulalization
from astropy.time import Time
from .aftmap import AFTmap
from numpy import ndarray, dtype

this_directory = Path(__file__).parent
__all__ = ['AFTmaps']


class AFTmaps:
    """
        A class for loading all AFT maps from a given directory or directories.


        Attributes
        ---------
        path: str | list[str]
            The path to the directory containing AFT map files.
        filetype: str
            The file extension of the AFT map files (e.g., "h5").
        date_fmt: str
            The date format string used to parse timestamps from filenames.
        filelist: list
             A list of file paths to AFT map files in the specified directory.
        filenames: numpy.ndarray
            An array of filenames extracted from the filelist.

        Parameters
        ----------
        path:
            The path or list of paths to the AFT map files. Defaults to the current working directory.
        filetype:
            The file extension of the AFT map. Default to the "aftmap"
        date_fmt:
            The date format string. Defaults to the "AFTmap_%Y%m%d_%H%M.h5"
        verbose:
            Whether to print or not. Defaults to the True
        hipft_prop:
            A dictionary for hipft properties, if filetype is "hipft", the hipft
    """

    def __init__(self, path: str | list | tuple = ".", filetype: str = "aftmap",
                 date_fmt: str = "AFTmap_%Y%m%d_%H%M.h5", monopole_corr: bool = False,
                 verbose: bool = True, hipft_prop: dict = None):
        """
        Initilaizing class function.
        """
        self.path = path
        self.filetype = filetype
        self.date_fmt = date_fmt
        self.verbose = verbose
        self.monopole_corr = monopole_corr
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

    def __repr__(self) -> str:
        return f"<AFTMaps(Path: {self.path}, Type: {self.filetype}, Total: {self.counts})>"

    def __str__(self) -> None:
        self.stats()

    @property
    def crn(self) -> np.ndarray:
        return self._CRN

    @property
    def aftmaps(self) -> np.ndarray:
        """
        An itterator containing all the AFTmaps associated.

        Yields
        ------
        aftmaps: AFTmaps
            AFTmaps object itterators.
        """
        for _file in self.filelist:
            yield AFTmap(_file, date_fmt=self.date_fmt, filetype=self.filetype)

    @property
    def timestamps(self) -> Time:
        """
        To get timestamps for all files in the AFT map directory.

        Returns
        -------
        timestamps: astropy.time.Time
            List of timestamps corresponding to each AFT map in the directory (astropy.time.Time).

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
        filelist: np.ndarray
            List of all files in the AFTmaps directory.
        """
        return self.filelist

    def get_filenames(self) -> np.ndarray:
        """
        Returns the names of the files in the AFTMap directory.

        Returns
        -------
        filebasename: np.ndarray
            List of files basenames.

        """
        return np.array(self.filenames)

    def stats(self) -> None:
        """
        Show the statistics of the loaded files.
        """
        _stats = {"RootDir": self.path, "FileType": self.filetype, "# Files": len(self.filelist),
                  "T-Initial": self.timestamps[0], "T-End": self.timestamps[-1]}
        for _s in _stats.keys():
            print("%-10s: %-30s" % (_s, _stats.get(_s)))

    def convert_all(self, convert_to: str = "fits", outpath: str = ".", **kwargs) -> None:
        """
        Convert all loaded AFT map files to the specified format.

        Parameters
        ----------
        convert_to: str, optional
            The output format to AFTMap to fits. Defaults to "fits"
        outpath: str, optional
            The directory in which to save the fits. Defaults to "."
        """
        for _file, i in zip(self.filelist, range(len(self.filelist))):
            mapobj = AFTmap(_file)
            yr, mo, day = mapobj.ymd
            _outpath = os.path.join(outpath, str(yr) + "/" + str(mo))
            _filename = os.path.splitext(os.path.basename(_file))[0] + "." + convert_to
            if not os.path.exists(_outpath):
                os.makedirs(_outpath, exist_ok=True)
            if convert_to == "png":
                mapobj.convert(convert_to=convert_to, verbose=self.verbose, outpath=_outpath, **kwargs)
            else:
                _filename = os.path.join(_outpath, _filename)
                mapobj.convert(convert_to=convert_to, verbose=self.verbose, outpath=_filename, **kwargs)
            if self.verbose:
                frac = float(i + 1) / len(self.filelist) * 100
                print(f"({frac:.2f}%) Converting {os.path.basename(_file)} to {os.path.basename(_filename)}", end="\r")

    @staticmethod
    def _generate_para(filelist: object) -> object:
        """
        A helper function to generate a aft parameters for a given file.
        Parameters
        ----------
        filelist: tuple
            The name of the file to generate the aft parameters.

        Returns
        -------
        para: tuple
            A tuple containing the all the AFT parameters for the file.

        """
        file, filetype, datefmt, monopole_corr = filelist
        mapobj = AFTmap(file, filetype=filetype, date_fmt=datefmt)
        tflux = round(np.abs(mapobj.flux).sum(), 3)
        pf = mapobj.polarfield(monopole_corr=monopole_corr)
        dp = mapobj.dipole(monopole_corr=monopole_corr)
        para = (mapobj.time, *pf, *dp, tflux)
        return para

    def generate_parameters(self, outfile: str = None, nthreds: int = None,
                            use_saved: bool = True) -> pd.DataFrame:
        """
        Function to generate AFT data for all the files.
        Parameters
        ----------
        use_saved: bool, optional
            Whether to relaod form saved file or not. Defaults to True.
        outfile: str, optional
            The file in which data is to be saved. Defaults to None.
        nthreds: int, optional
            The number of thredds to use for the calculations. Defaults to None.

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
                if outfile is not None:
                    df.to_csv(outfile, float_format="%.3g", index=False)
                return df
            elif self.verbose:
                print("WARNING: Saved data is outdated, skipping saved data.")

        result = []
        tint = time.time()
        bar_format = '{desc}: {percentage:3.2f}%|{bar}{r_bar}'
        inputpara = [(file, self.filetype, self.date_fmt, self.monopole_corr) for file in self.filelist]
        with Pool(nthreds) as pool:
            for _para in tq.tqdm(pool.imap(self._generate_para, inputpara),
                                 bar_format=bar_format, total=self.counts):
                result.append(_para)
        if self.verbose:
            print(f"Completed in {int(time.time() - tint)} seconds.")
            print(f"AFT parameters are saved to {outfile}.")
        df = pd.DataFrame(result, columns=["Time", "PolarN", "PolarS", "ADM", "EDM", "TotalFlux"])
        if outfile is not None:
            df.to_csv(outfile, float_format="%.3g", index=False)
        df2h5(df, save_file)
        return df

    def get_crmap(self, cr_number: int) -> np.ndarray:
        """

        Parameters
        ----------
        cr_number:
            Carrington Map for a given Carrington Rotation Number.

        Returns
        -------
        crmap: np.ndarray
            Carrington Map for a given Carrington Rotation Number.

        """
        if cr_number not in self._CRN:
            raise ValueError("The given Carrington Rotation is Not present in the data.")

        ind = np.where(self._crn == cr_number)
        crmap = 0
        for _fl in self.filelist[ind]:
            aftobj = AFTmap(_fl, filetype=self.filetype, date_fmt=self.date_fmt)
            crmap = crmap + aftobj.aftmap
        crmap = crmap / ind[0].size
        return crmap

    def cravgmap(self, use_saved: bool = True) -> np.ndarray:
        """

        Parameters
        ----------
        use_saved:
            Whether to relaod form saved file or not. Defaults to True.

        Returns
        -------
        crmap: np.ndarray
            Carrington Rotation Averaged Map for all the maps. Basically Butterfly Diagram.

        """
        save_file = os.path.join(this_directory, ".bflydata.h5")
        save_exist = os.path.isfile(save_file)
        if use_saved & save_exist:
            with hdf.File(save_file, "r") as fh5:
                bflymap = fh5["bdata"][()]
            if bflymap.shape[1] == len(self._CRN):
                return bflymap
            elif self.verbose:
                print("WARNING: Saved data is outdated, skipping saved data.")

        # Get carrington Numbers
        crn = self._CRN
        crn0 = crn.min()
        bar_format = '{desc}: {percentage:3.2f}%|{bar}{r_bar}'
        bflymap = np.zeros((512, len(crn)))
        for _crn in tq.tqdm(crn, bar_format=bar_format, total=len(crn)):
            bflymap[:, _crn - crn0] = self.get_crmap(_crn).mean(axis=1)
        if save_exist:
            os.remove(save_file)
        with hdf.File(save_file, "w") as fh5:
            fh5.create_dataset("bdata", data=bflymap)

        return bflymap

    @property
    def visualize(self) -> Visulalization:
        """

        Returns
        -------
        Visulalization: Visulalization
            AFTpy Visulalization object.

        """
        return Visulalization(self)
