"""
SFT Module: A Python Library for Handling SFT Magnetic Field Data

This module provides classes and methods for managing and analyzing
magnetic field data from the Advanced Flux Telescope (AFT). The SFT
class initializes with a list of relevant AFT map files and provides
methods for processing them, generating AFT parameters, and creating
butterfly diagrams.

Classes
-------
- JsonSerialize : Custom JSON encoder for handling NumPy data types.
- SFT : Handles AFT map data, iterates over SFTMap objects, and generates
        AFT parameters and butterfly diagrams.

Features
--------
- Loads AFT map data from single files, directories, or direct data-header pairs.
- Iterates over AFT map files as `SFTMap` objects.
- Computes and stores solar magnetic parameters such as:
  - Total absolute flux
  - Axial and equatorial dipole moments
  - Polar region magnetic strengths
  - Azimuthal averages
- Supports parallel processing for database generation.
- Exports results in JSON, ASDF, or HDF5 formats.

Usage Example
-------------
>>> from sft import SFT
>>> sft = SFT("/path/to/data")
>>> for aft_map in sft:
>>>     print(aft_map.parameters())

References
----------
* AFT Mission Paper
* SunPy Library: https://sunpy.org/
"""

import os
import json
import pathlib
import datetime as dt
from pathlib import Path
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool as Pool
from typing import Any

import numpy as np
from asdf import AsdfFile
from tqdm.autonotebook import tqdm

from .utilities import file_search
from .sftbasemap import SFTMap
from .io.util import is_path, is_url, is_dir

this_directory = Path.home() / ".aftpy"
__all__ = ["SFT"]


class JsonSerialize(json.JSONEncoder):
    """
    Custom JSON encoder for handling NumPy types.

    Methods
    -------
    default(obj)
        Converts NumPy types to Python native types for JSON serialization.
    """

    def default(self, obj):
        """
        Convert NumPy data types to native Python types.

        Parameters
        ----------
        obj : object
            The object to be serialized.

        Returns
        -------
        object
            Converted object suitable for JSON encoding.
        """
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, str):
            return obj.strip()
        return super().default(obj)


class SFT:
    """
    Handles AFT map data by initializing with a list of files and
    providing methods to process them.

    Parameters
    ----------
    *args : tuple
        Can be:
        - A single file path or URL
        - A directory containing multiple files
        - A data-header pair for direct SFTMap creation
    **kwargs : dict
        Additional parameters for custom configurations.

    Attributes
    ----------
    nfiles : int
        Number of AFT map files detected.
    files : list of str
        List of file paths.
    sftmodel : str or None
        Detected SFT model type.
    """

    def __init__(self, *args, **kwargs):
        self.nargs = len(args)
        self.kwargs = kwargs
        self.nfiles = 0
        self.current = -1
        self.sftmodel = None
        self.files = None

        if self.nargs == 1:
            self.input = args[0]

            # Single file case (Path or URL)
            if isinstance(self.input, str) and (is_path(self.input) or is_url(self.input)):
                self.nfiles = 1
                self.sftmodel = "aftv2"

            # Directory case
            elif isinstance(self.input, str) and is_dir(self.input):
                file_patterns = {
                    "aftv2": "AFTmap*.h5",
                    "others": "*.fits",
                    "aftv1": "*.dat"
                }

                for model, pattern in file_patterns.items():
                    files = file_search(self.input, fileext=pattern)
                    if len(files) > 0:
                        self.sftmodel = model
                        self.files = files
                        self.nfiles = len(files)
                        break
                    else:
                        raise ValueError(f"No AFTmap or SFTmap of known types found in {self.input}.")

            else:
                raise ValueError(f"Invalid input: {self.input}")

        elif self.nargs == 2:
            self._data = args[0]
            self._header = args[1]
            self.nfiles = 1
        else:
            raise ValueError("Unsupported number of arguments.")

    def __len__(self):
        return self.nfiles

    def __iter__(self):
        """Returns an iterator instance."""
        return self

    def __next__(self):
        """
        Returns the next SFTMap object.

        Returns
        -------
        SFTMap
            A processed SFTMap instance.

        Raises
        ------
        StopIteration
            If there are no more files to process.
        """
        self.current += 1
        if self.current >= self.nfiles:
            raise StopIteration("Max number of files reached.")

        if self.nfiles == 1:
            if self.nargs == 2:
                return SFTMap(self._data, self._header)
            elif self.nargs == 1:
                return SFTMap(self.input)

        elif self.nfiles > 1 and self.nargs == 1:
            return SFTMap(self.files[self.current])

        # return None

    @staticmethod
    def __get_para(x):
        return x.parameters()

    def generate_database(self, azavg=True, outfile=None, n_threads=None, parallel=True, dbfmt="JSON"):
        """
        Generates a database of AFT parameters for all files.

        Parameters
        ----------
        azavg : bool, optional
            Whether to include azimuthal averages.
        outfile : str,Path, None, optional
            Output file path.
        n_threads : int, optional
            Number of threads for parallel processing.
        parallel : bool, optional
            Whether to use parallel processing.
        dbfmt : str, optional
            Database format ("JSON", "ASDF", or "HDF5").

        Returns
        -------
        None, dict
            Saves the results in the specified format.
        """
        if n_threads is None:
            n_threads = max(1, cpu_count() - 1)

        ext_map = {"JSON": ".json", "ASDF": ".asdf", "HDF5": ".hdf5"}
        results: dict[Any, Any] = {}
        if dbfmt.upper() != "DICT":
            if outfile is None:
                outfile = f"SFT_db_{dt.datetime.now().strftime('%Y%m%d')}{ext_map[dbfmt.upper()]}"
                path = Path(outfile)
            elif isinstance(outfile, Path):
                path = outfile
            elif isinstance(outfile, str):
                path = Path(outfile)
            else:
                raise ValueError(f"Invalid output file type: {type(outfile)}")

            if path.is_dir():
                outfile = os.path.join(outfile, f"SFT_db_{dt.datetime.now().strftime('%Y%m%d')}{ext_map[dbfmt.upper()]}")
            else:
                dbfmt = path.suffix[1:].upper()

        if parallel:
            with Pool(n_threads) as pool:
                for result in tqdm(pool.imap(self.__get_para, self),
                                   total=len(self.files)-self.current-1):
                    if not azavg:
                        result.pop("AZAVG")
                    results[result.get("DATE_OBS", None)] = result
        else:
            for aft_map in self:
                result = self.__get_para(aft_map)
                results[result.get("DATE_OBS", None)] = result

        # Save database
        if dbfmt.upper() == "JSON":
            with open(outfile, "w") as f:
                json.dump(results, f, indent=4, cls=JsonSerialize)
            print(f"Wrote {len(results)} files to {outfile}.")
            return None
        elif dbfmt.upper() == "ASDF":
            AsdfFile(results).write_to(outfile)
            print(f"Wrote {len(results)} files to {outfile}.")
            return None
        elif dbfmt.upper() == "DICT":
            print(f"Wrote {len(results)} files to {outfile}.")
            return results
        else:
            print(f"Unrecognized database format: {dbfmt} returning dictionary",)
            return results


    def upadate_database(self, dbfile, verbose=True):
        """

        Parameters
        ----------
        verbose
        dbfile

        Returns
        -------

        """
        outfile = Path(dbfile)
        tempfile = Path(".temp.json")
        if not outfile.exists():
            FileNotFoundError(f"{outfile} does not exist.")
        with open("/Users/bjha/Data/AFT/AFTmapv2/SFT_db_20250320.json") as fl:
            results = json.load(fl)
        if verbose:
            print(f"Updating database of {self.nfiles-len(results)} files.")
        if len(results) == self.nfiles:
            print(f"Based on {len(results)} files databse seems updated.")
            return None
        self.current = len(results)-1
        if verbose:
            print(f"Current file: {self.current} / {self.nfiles} files.")
        _results = self.generate_database(dbfmt="dict")
        if _results:
            results.update(_results)
        with open(tempfile, "w") as f:
            json.dump(results, f, indent=4, cls=JsonSerialize)
        outfile.unlink()
        tempfile.rename(outfile)
        print(f"Updated {outfile} files.")
        return None
