import os
import numpy as np
import h5py as h5
import pandas as pd


def file_search(path: str, fileext: str = None) -> np.ndarray:
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
            elif f.endswith(fileext) and (not f.startswith(".")):
                files.append(os.path.join(dr, f))
    return np.sort(np.array(files))


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
