# AFTpy
The Python package enables users to download, visualize, and analyze
AFT HDF (h5) data seamlessly. `AFTPy` offers a comprehensive solution for
analyzing and downloading Advective Flux Transport (AFT) data,
streamlining the process of data downloading in parallel and
also facilitating the conversion of H5 files into the most popular
FITS file format or various other formats with ease.

## Installation
### From PyPI
```shell
pip install aftpy
```

### From source
```shell
git clone git@github.com:bibhuraushan/aftpy.git
cd aftpy
python setup.py install
```

## Discriptions
`aftpy` provides two importent python module named `aftmap` and `aftgetdata`. These modules can be used
to read the single **aftmap file** or load all the files from a given directory.

## aftmap module
`aftmap` modules also provides two python class `AFTmap` and `AFTload`. The `AFTmap` class
provide an interface to read the singale H5 **AFTmap** file and provide you the functions and isntances 
to get the information and plot the data. The other class `AFTload` provides the interface to load 
all the data from a directory and provide the instances and function to know about the loaded data. It also
provides a function to convert all the loaded data in to popular **FITS** file.

### AFTmap Class

### AFTload Class

A class for loading all AFT maps from directory.

#### Attributes

- `path` (str): The path to the directory containing AFT map files.
- `filetype` (str): The file extension of the AFT map files (e.g., "h5").
- `date_fmt` (str): The date format string used to parse timestamps from filenames.
- `filelist` (list): A list of file paths to AFT map files in the specified directory.
- `filenames` (numpy.ndarray): An array of filenames extracted from the filelist.

#### Methods

- `convert_all(convert_to="fits", outpath=".", verbose=True)` Convert all loaded AFT map files to the specified format.
  - `convert_to` (str, optional): The output format to convert the AFT map files to. Defaults to "fits".
  - `outpath` (str, optional): The directory path to save the converted files. Defaults to current directory.
  - `verbose` (bool, optional): Whether to print conversion progress. Defaults to True.

#### Example Usage

```python
import aftpy.aftmap as aft

# Initialize AFTload object
loader = aft.AFTload(path="/path/to/aft/maps", filetype="h5")

# Convert all AFT map files to FITS format to '/path/to/converted/maps'
loader.convert_all(convert_to="fits", outpath="/path/to/converted/maps", verbose=True)
```


## aftgetdata module

### AFTdownlaod Clas
A class for downloading AFT map files from a specified URL.

#### Attributes

- `ncpu` (int): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.
- `dt` (module): Alias for the `datetime` module.
- `root_url` (str): The root URL from where AFT map files will be downloaded. Defaults to "https://data.boulder.swri.edu/lisa/".
- `urls` (list): List of URLs of AFT map files.
- `datafile` (str): File path to store the list of files in CSV format.
- `datalist` (DataFrame): DataFrame containing the list of files and corresponding timestamps.

#### Methods
- `get_list(t0=None, t1=None, dt=1) -> data (DataFrame)`
  - `t0` (datetime.datetime, optional): Start time of the time range. Defaults to None.
  - `t1` (datetime.datetime, optional): End time of the time range. Defaults to None.
  - `dt` (int, optional): Time interval for sampling files within the time range. Defaults to 1.
  - `Returns` data (DataFrame): DataFrame containing the list of files within the specified time range.

- `reload_files(url=None, filetype="h5")` Reload the list of AFT map files from the root URL.
  - `url` (str, optional): The URL to reload the list of files from. Defaults to None.
  - `filetype` (str, optional): The file extension of AFT map files. Defaults to "h5".
  - `Returns`: True if the list of files is successfully reloaded.
- `download(dataframe, rootpath=None, ncpu=None)` Download AFT map files listed in the DataFrame.
  - `dataframe` (DataFrame): DataFrame containing the list of files to download.
  - `rootpath` (str, optional): Root directory path to save the downloaded files. Defaults to None.
  - `ncpu` (int, optional): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.


#### Example Usage

```python
import aftpy.getaftdata as aftget
# Initialize AFTdownload object
downloader = aftget.AFTdownload()

# Reload the list of AFT map files
downloader.reload_files()

# Get the list of AFT map files within a specified time range
file_list = downloader.get_list(t0=dt.datetime(2023, 1, 1), t1=dt.datetime(2023, 1, 7))

# Download AFT map files listed in the DataFrame
downloader.download(file_list)
```
