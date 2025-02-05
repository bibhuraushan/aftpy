"""

"""
import h5py as h5
from astropy.io import fits
from astropy.io.fits.header import Card
from .util import is_url, is_path

def read_aftv2(filename: str) -> dict:
    """
    Reads AFT v2 data from an HDF5 file and stores it in an SFTdata object.

    Parameters
    ----------
    filename : str
        Path to the HDF5 file containing AFT v2 data.

    Returns
    -------
    SFTdata
        An object containing the AFT map data, mask, and FITS header.
    """
    sft_data = {}
    if is_path(filename):
        _filename = filename
    elif is_url(filename):
        _filename = filename
    else:
        raise FileNotFoundError(f"File not found: {filename}")

    with h5.File(filename, "r") as fl:
        sft_data["data"] = fl["/maps/aftmap"][()]
        try:
            sft_data["mask"] = fl["/maps/mask"][()]
        except KeyError:
            sft_data["mask"] = None
        try:
            sft_data["magmap"] = fl["/maps/magmap"][()]
        except KeyError:
            sft_data["magmap"] = None
        header = []
        for _key in fl["header"].keys():
            _val = list(fl["header"][_key])[0]
            if isinstance(_val, bytes):
                _val = _val.decode("utf-8")
            try:
                _info = list(fl["header"][_key].attrs["descr"])[0].decode("utf-8")
            except KeyError:
                _info = "None"
            header.append(Card(_key, _val, _info))

        sft_data["header"] = fits.Header(header)

    return sft_data