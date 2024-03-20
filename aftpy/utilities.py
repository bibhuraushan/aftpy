import os
import numpy as np


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
            elif f.endswith(fileext):
                files.append(os.path.join(dr, f))
    return np.sort(np.array(files))
