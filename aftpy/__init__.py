"""
The `aftmap` module also provides two Python classes `AFTmap` and `AFTmaps` (earlier `AFTload`). The `AFTmap` class
provides an interface to read a single H5 **AFTmap** file and provides you the functions and instances
to get the information and plot the data. The other class `AFTmaps` provides the interface to load
all the data from a directory and provide the instances and function to know about the loaded data. It also
provides a function to convert all the loaded data into **FITS** as well as many other formats.

"""
import importlib.metadata
__version__ = importlib.metadata.version("aftpy")

__author__ = 'Bibhuti Kumar Jha'
__email__ = 'bibhuraushan1@gmail.com'
__all__ = ['AFTmap', 'AFTmaps', 'AFTdownload']

from .aftmap import AFTmap
from .aftmaps import AFTmaps
from .getaftdata import AFTdownload
from .visulalization import Visulalization


