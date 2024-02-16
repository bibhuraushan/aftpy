import datetime as dt
import requests
from urllib.request import Request, urlopen, urlretrieve
from bs4 import BeautifulSoup
import pandas as pd
import os
import numpy as np
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
import tqdm as tq
from pathlib import Path

this_directory = Path(__file__).parent


def _download(args):
    """
    Internal function to download AFT map files.

    Args:
        args (tuple): Tuple containing URL and root path of the file to download.

    Returns:
        str: Message indicating whether the file was successfully downloaded or not.
    """
    url, path = args[0], args[1]
    filename = os.path.basename(url)
    yr, mo = filename[7:11], filename[11:13]
    fullpath = os.path.join(path, f"{str(yr)}/{str(mo)}")
    if not os.path.exists(fullpath):
        os.makedirs(fullpath, exist_ok=True)
    fullpath = os.path.join(fullpath, filename)
    try:
        urlretrieve(url, fullpath)
        return f"{os.path.basename(url)} sucessfully downloded."
    except:
        print(f"Unable to download file {os.path.basename(url)}")


class AFTdownload:
    """
       A class for downloading AFT map files from a specified URL.

       Attributes:
           ncpu (int): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.
           dt (module): Alias for the `datetime` module.
           root_url (str): The root URL from where AFT map files will be downloaded.
           urls (list): List of URLs of AFT map files.
           times (list): List of timestamps corresponding to the AFT map files.
           counts (int): Total count of AFT map files.
           datafile (str): File path to store the list of files in CSV format.
           datalist (DataFrame): DataFrame containing the list of files and corresponding timestamps.

       Methods:
           __init__(root_url="https://data.boulder.swri.edu/lisa/"):
               Initializes the AFTdownload object.

           get_list(t0=None, t1=None, dt=1):
               Gets a list of AFT map files within a specified time range.

           reload_files(url=None, filetype="h5"):
               Reloads the list of AFT map files from the root URL.

           download(dataframe, rootpath=None, ncpu=None):
               Downloads AFT map files listed in the DataFrame.

           _download(args):
               Internal function to download AFT map files.
       """
    ncpu = cpu_count() - 1

    def __init__(self, root_url: str = "https://data.boulder.swri.edu/lisa/"):
        """
        Initializes the AFTdownload object.

        Args:
            root_url (str, optional): The root URL from where AFT map files will be downloaded.
                Defaults to "https://data.boulder.swri.edu/lisa/".
        """
        self.root_url = root_url
        self.urls = []
        self.times = []
        self.counts = 0
        self.datafile = f"{this_directory}/list_of_files.csv"
        if os.path.exists(self.datafile):
            self.datalist = pd.read_csv(self.datafile, index_col="times")
        else:
            self.reload_files()
            self.datalist = pd.read_csv(self.datafile, index_col="times")

    def get_list(self, t0: dt.datetime = None, t1: dt.datetime = None, cadance: int = 1) -> pd.DataFrame:
        """
                Gets a list of AFT map files within a specified time range.

                Args:
                    t0 (datetime.datetime, optional): Start time of the time range. Defaults to None.
                    t1 (datetime.datetime, optional): End time of the time range. Defaults to None.
                    dt (int, optional): Time interval for sampling files within the time range. Defaults to 1.

                Returns:
                    DataFrame: DataFrame containing the list of files within the specified time range.
                    :param cadance:
                    :param t0:
                    :param t1:
                """
        deltat = int(4 / cadance)
        if (t0 is None) & (t1 is None):
            return self.datalist
        else:
            t0, t1 = t0.isoformat(), t1.isoformat()
            ind0 = np.where(self.datalist.index == t0)
            ind1 = np.where(self.datalist.index == t1)
            if ind0[0].size != 0:
                x1 = ind0[0][0]
            else:
                x1 = 0

            if ind1[0].size != 0:
                x2 = ind1[0][0]
            else:
                x2 = -1
            data = self.datalist.iloc[x1:x2:deltat].copy()
            return data

    def reload_files(self, url: str = None, filetype: str = "h5") -> bool:
        """
        Reloads the list of AFT map files from the root URL.

        Args:
            url (str, optional): The URL to reload the list of files from. Defaults to None.
            filetype (str, optional): The file extension of AFT map files. Defaults to "h5".

        Returns:
            bool: True if the list of files is successfully reloaded.
        """
        if url is None:
            url = self.root_url
        url = url.replace(" ", "%20")
        req = Request(url)
        a = urlopen(req).read()
        soup = BeautifulSoup(a, 'html.parser')
        x = (soup.find_all('a'))
        for i in x:
            file_name = i.extract().get_text()
            url_new = url + file_name
            url_new = url_new.replace(" ", "%20")
            if file_name.endswith(filetype):
                self.urls.append(url_new)
                _time = dt.datetime.strptime(file_name, "AFTmap_%Y%m%d_%H%M.h5").isoformat()
                self.times.append(_time)
                self.counts += 1
                print(f"Total {self.counts} files located.", end="\r")
            if file_name[-1] == '/' and file_name[0] != '.':
                self.reload_files(url=url_new)
        data = pd.DataFrame({"times": self.times, "urls": self.urls, })
        data.set_index("times", inplace=True)
        data.to_csv(self.datafile)
        return True

    def download(self, dataframe, rootpath: str = None, ncpu: int = None):
        """
        Downloads AFT map files listed in the DataFrame.

        Args:
            dataframe (DataFrame): DataFrame containing the list of files to download.
            rootpath (str, optional): Root directory path to save the downloaded files. Defaults to None.
            ncpu (int, optional): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.
        """
        if ncpu is None:
            self.ncpu = cpu_count() - 1
        else:
            self.ncpu = ncpu
        if rootpath is None:
            os.mkdir("data")
            rootpath = "data"
        nfiles = dataframe.shape[0]
        paths = np.full(nfiles, rootpath)
        urls = dataframe.urls.values
        args = zip(urls, paths)
        pool = ThreadPool(self.ncpu)
        for result in tq.tqdm(pool.imap_unordered(_download, args), total=urls.size):
            pass
        pool.close()
        pool.join()
