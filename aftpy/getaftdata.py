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
    ncpu = cpu_count() - 1

    def __init__(self, root_url="https://data.boulder.swri.edu/lisa/"):
        self.dt = dt
        self.root_url = root_url
        self.urls = []
        self.times = []
        self.counts = 0
        self.datafile = f"{this_directory}/list_of_files.csv"
        if os.path.exists(self.datafile):
            self.datalist = pd.read_csv(self.datafile, index_col="times")
        else:
            self.reload_files()

    def get_list(self, t0=None, t1=None, dt=1):
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
            data = self.datalist.iloc[x1:x2:dt].copy()
            return data

    def reload_files(self, url=None, filetype="h5"):
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

    def download(self, dataframe, rootpath=None, ncpu=None):
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
