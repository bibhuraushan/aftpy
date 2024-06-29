![aftlogo](aftpy_logo.png)
# AFTpy

`AFTPy` offers a comprehensive solution for
analyzing and downloading Advective Flux Transport (AFT) data,
streamlining the process of data downloading in parallel and
also facilitating the conversion of H5 files into the most popular
FITS file format or various other formats with ease.


---

`aftpy` provides two important Python modules named `aftmap` and `aftgetdata`. These modules can be used
to read a single **aftmap file** or load all the files from a given directory.

## aftmap module

---

The `aftmap` module also provides two Python classes `AFTmap` and `AFTmaps` (earlier `AFTload`). The `AFTmap` class
provides an interface to read a single H5 **AFTmap** file and provides you the functions and instances 
to get the information and plot the data. The other class `AFTmaps` provides the interface to load 
all the data from a directory and provide the instances and function to know about the loaded data. It also
provides a function to convert all the loaded data into **FITS** as well as many other formats.

## aftgetdata module

___

The `aftgetdata` module provide a one step approach to download `AFTmap`
file in HDF file format from the  SwRI Public Data Server. This uses a parralel data download
approach to download data very quickely. 
