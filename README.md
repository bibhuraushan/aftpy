# AFTpy
An installable python package to download, visualize and analyze AFT HDF (h5) data.

## aftmap
This is the module that read
either one file or all files from a given directory as an object.

These objects can be used to visulaize the data use tha data
and also convert the data to desired format.

## AFTload Class

A class for loading AFT maps from files.

### Attributes

- `path` (str): The path to the directory containing AFT map files.
- `filetype` (str): The file extension of the AFT map files (e.g., "h5").
- `date_fmt` (str): The date format string used to parse timestamps from filenames.
- `filelist` (list): A list of file paths to AFT map files in the specified directory.
- `filenames` (numpy.ndarray): An array of filenames extracted from the filelist.

### Methods

### `__init__(path=".", filetype="h5", date_fmt="AFTmap_%Y%m%d_%H%M.h5", verbose=True)`

Initialize the AFTload object.

#### Arguments

- `path` (str, optional): The path to the directory containing AFT map files. Defaults to current directory.
- `filetype` (str, optional): The file extension of the AFT map files. Defaults to "h5".
- `date_fmt` (str, optional): The date format string used to parse timestamps from filenames. Defaults to "AFTmap_%Y%m%d_%H%M.h5".
- `verbose` (bool, optional): Whether to print statistics about the loaded AFT map files. Defaults to True.

### `convert_all(convert_to="fits", outpath=".", verbose=True)`

Convert all loaded AFT map files to the specified format.

#### Arguments

- `convert_to` (str, optional): The output format to convert the AFT map files to. Defaults to "fits".
- `outpath` (str, optional): The directory path to save the converted files. Defaults to current directory.
- `verbose` (bool, optional): Whether to print conversion progress. Defaults to True.

### Example Usage

```python
# Initialize AFTload object
loader = AFTload(path="/path/to/aft/maps", filetype="h5")

# Convert all AFT map files to FITS format
loader.convert_all(convert_to="fits", outpath="/path/to/converted/maps", verbose=True)
```

# AFTdownload Class

A class for downloading AFT map files from a specified URL.

## Attributes

- `ncpu` (int): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.
- `dt` (module): Alias for the `datetime` module.
- `root_url` (str): The root URL from where AFT map files will be downloaded.
- `urls` (list): List of URLs of AFT map files.
- `times` (list): List of timestamps corresponding to the AFT map files.
- `counts` (int): Total count of AFT map files.
- `datafile` (str): File path to store the list of files in CSV format.
- `datalist` (DataFrame): DataFrame containing the list of files and corresponding timestamps.

## Methods

### `__init__(root_url="https://data.boulder.swri.edu/lisa/")`

Initialize the AFTdownload object.

#### Arguments

- `root_url` (str, optional): The root URL from where AFT map files will be downloaded. Defaults to "https://data.boulder.swri.edu/lisa/".

### `get_list(t0=None, t1=None, dt=1)`

Get a list of AFT map files within a specified time range.

#### Arguments

- `t0` (datetime.datetime, optional): Start time of the time range. Defaults to None.
- `t1` (datetime.datetime, optional): End time of the time range. Defaults to None.
- `dt` (int, optional): Time interval for sampling files within the time range. Defaults to 1.

#### Returns

- `data` (DataFrame): DataFrame containing the list of files within the specified time range.

### `reload_files(url=None, filetype="h5")`

Reload the list of AFT map files from the root URL.

#### Arguments

- `url` (str, optional): The URL to reload the list of files from. Defaults to None.
- `filetype` (str, optional): The file extension of AFT map files. Defaults to "h5".

#### Returns

- `bool`: True if the list of files is successfully reloaded.

### `download(dataframe, rootpath=None, ncpu=None)`

Download AFT map files listed in the DataFrame.

#### Arguments

- `dataframe` (DataFrame): DataFrame containing the list of files to download.
- `rootpath` (str, optional): Root directory path to save the downloaded files. Defaults to None.
- `ncpu` (int, optional): Number of CPU cores to utilize for downloading files. Defaults to `cpu_count() - 1`.

### `_download(args)`

Internal function to download AFT map files.

#### Arguments

- `args` (tuple): Tuple containing URL and root path of the file to download.

#### Returns

- `str`: Message indicating whether the file was successfully downloaded or not.

## Example Usage

```python
# Initialize AFTdownload object
downloader = AFTdownload()

# Reload the list of AFT map files
downloader.reload_files()

# Get the list of AFT map files within a specified time range
file_list = downloader.get_list(t0=dt.datetime(2023, 1, 1), t1=dt.datetime(2023, 1, 7))

# Download AFT map files listed in the DataFrame
downloader.download(file_list)
```
