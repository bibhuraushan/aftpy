import os
import numpy as np
import h5py as h5
import pandas as pd
import re
import os
from datetime import datetime
from fnmatch import fnmatch
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all_ = ["datetime_from_filename", "file_search"]

def datetime_from_filename(filename, isoformat=False):
    """
    Extracts date and time information from a filename and returns a datetime object.

    The function scans the filename for common date-time patterns and attempts to
    extract the date-time in various formats, such as:
    - YYYY-MM-DD_HH:MM:SS
    - YYYYMMDD_HHMM
    - YYMMDD_HHMMSS
    - HH:MM:SS
    - DD-MM-YYYY_HH:MM:SS

    If no date-time information is found, the function returns `None`.

    Parameters
    ----------
    filename : str
        The input filename string, typically containing date-time information.
    isoformat : bool, optional
        If True, returns the extracted datetime as an ISO 8601 formatted string (default is False).

    Returns
    -------
    datetime or str or None
        - `datetime.datetime` object if extraction is successful.
        - `str` in ISO 8601 format (`YYYY-MM-DDTHH:MM:SS`) if `isoformat=True`.
        - `None` if no valid date/time is found in the filename.
    """

    # Define regex patterns with corresponding datetime formats
    patterns = [
        (r'(\d{4})[-_]?(\d{2})[-_]?(\d{2})[T_]?(\d{2})[:_]?(\d{2})[:_]?(\d{2})', "%Y%m%d%H%M%S"),  # YYYY-MM-DDTHHMMSS
        (r'(\d{4})[-_]?(\d{2})[-_]?(\d{2})[T_]?(\d{2})[:_]?(\d{2})', "%Y%m%d%H%M"),  # YYYY-MM-DD_HH:MM
        (r'(\d{4})[-_]?(\d{2})[-_]?(\d{2})', "%Y%m%d"),  # YYYY-MM-DD
        (r'(\d{8})[T_]?(\d{6})', "%Y%m%d%H%M%S"),  # YYYYMMDD_HHMMSS
        (r'(\d{8})[T_]?(\d{4})', "%Y%m%d%H%M"),  # YYYYMMDD_HHMM
        (r'(\d{8})', "%Y%m%d"),  # YYYYMMDD
        (r'(\d{6})[T_]?(\d{6})', "%y%m%d%H%M%S"),  # YYMMDD_HHMMSS
        (r'(\d{2})[-_](\d{2})[-_](\d{4})[T_](\d{2}):(\d{2}):(\d{2})', "%d%m%Y%H%M%S"),  # DD-MM-YYYY_HH:MM:SS
        (r'(\d{6})[T_]?(\d{4})', "%y%m%d%H%M"),  # YYMMDD_HHMM
        (r'(\d{6})', "%y%m%d"),  # YYMMDD
        (r'(\d{2}):(\d{2}):(\d{2})', "%H:%M:%S"),  # HH:MM:SS
        (r'(\d{2}):(\d{2})', "%H:%M"),  # HH:MM
    ]

    _filename = os.path.basename(filename)  # Ensure only filename is processed

    for pattern, date_format in patterns:
        match = re.search(pattern, _filename)
        if match:
            try:
                extracted_time = "".join(match.groups())  # Join matched parts without extra spaces
                dt = datetime.strptime(extracted_time, date_format)
                return dt.isoformat() if isoformat else dt  # Return ISO format if requested
            except ValueError:
                continue  # Skip if parsing fails

    return None  # Return None if no date-time is found



def file_search(path: str, fileext: str = None) -> np.ndarray|None:
    """
       Search for files with the specified file extension in the given directory.

    Parameters
    ----------
    path: str
    The directory in which the given file will be searched.
    fileext: str
    The file extension to be searched for. Defaults to "" (no filtering)
    """
    files = []
    for dr, _, file in os.walk(path):
        for f in file:
            if fileext is None:
                files.append(os.path.join(dr, f))
            elif fnmatch(f, fileext):
                files.append(os.path.join(dr, f))
    if len(files) != 0:
        return np.sort(np.array(files))
    else:
        return None

def detect_sftmodel(string: str) -> str | None:
    """

    Parameters
    ----------
    string

    Returns
    -------

    """
    if os.path.basename(string).lower().startswith("aftmap"):
        sft_model = "aftv2"
        logger.debug(f"AFTmap detected: {sft_model}")
        return sft_model
    elif os.path.basename(string).lower().startswith("adapt"):
        sft_model = "adapt"
        logger.debug(f"Adapt detected: {sft_model}")
        return sft_model
    else:
        sft_model = None
        logger.debug(f"Unknown sftmodel detected: {sft_model}")
        return sft_model


def df2h5(df: pd.DataFrame, filename: str) -> None:
    if os.path.isfile(filename):
        os.remove(filename)
    with h5.File(filename, "w") as fl:
        fl.create_dataset("time", data=df.Time.values.astype("S"),)
        fl.create_dataset("pfn", data=df.PolarN.values)
        fl.create_dataset("pfs", data=df.PolarS.values)
        fl.create_dataset("adm", data=df.ADM.values)
        fl.create_dataset("edm", data=df.EDM.values)
        fl.create_dataset("tflux", data=df.TotalFlux.values)


def h52df(filename: str) -> pd.DataFrame:
    if os.path.isfile(filename):
        df = pd.DataFrame({})
        with h5.File(filename, "r") as fl:
            df["Time"] = np.char.decode(fl["time"][()], "UTF-8")
            df["PolarN"] = fl["pfn"][()]
            df["PolarS"] = fl["pfs"][()]
            df["ADM"] = fl["adm"][()]
            df["EDM"] = fl["edm"][()]
            df["TotalFlux"] = fl["tflux"][()]
    else:
        raise FileNotFoundError(f"File {filename}")

    return df
