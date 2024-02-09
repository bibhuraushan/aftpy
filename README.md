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
