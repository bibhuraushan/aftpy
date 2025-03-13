import numpy as np
import pylab as p
from matplotlib import pyplot as plt
from matplotlib.axis import Axis

import numpy as np
import matplotlib.pyplot as plt
from .util import cmap_balance, cmap_YlGrBl
from pathlib import Path
import asdf
import pandas as pd


def plot_butterfly(data, ftime, crn=None, sinlat=True,
                   ax=None, plot_prop=None, colorbar_prop=None):
    """
    Plots the magnetic butterfly diagram, representing the latitudinal distribution of magnetic fields over time.
    Optionally, a secondary axis can be added to display Carrington rotation numbers.

    Parameters
    ----------
    data : ndarray
        2D array representing the magnetic field data (latitude vs. time).
    ftime : ndarray
        1D array of time values corresponding to the x-axis.
    crn : ndarray, optional
        1D array of Carrington rotation numbers (same length as `time`).
    sinlat : bool, optional, default=True
        If True, plots against sine(latitude); otherwise, uses latitude directly.
    ax : matplotlib.axes._subplots.Axes, optional
        Axes object to plot on. If None, a new figure is created.
    plot_prop : dict, optional
        Dictionary of properties passed to `ax.pcolormesh()` for customizing the plot.
    colorbar_prop : dict, optional
        Dictionary of properties passed to `fig.colorbar()` for customizing the colorbar.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The primary axes object of the plot.

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> time = np.linspace(2000, 2025, 100)
    >>> latitudes = np.linspace(-90, 90, 50)
    >>> data = np.random.uniform(-10, 10, (50, 100))
    >>> crn = np.linspace(1950, 2250, 100)  # Simulated Carrington Rotation Numbers
    >>> fig, ax = plot_butterfly(data, ftime, crn)
    >>> plt.show()
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure  # Get figure from the provided axes

    # Define latitude range
    latitudes = np.linspace(-90, 90, data.shape[0])
    if sinlat:
        y_vals = np.sin(np.radians(latitudes))
        y_label = "Sine Latitude"
    else:
        y_vals = latitudes
        y_label = "Latitude (°)"

    # Default plotting properties
    plot_kwargs = {"cmap": cmap_YlGrBl(), "vmax": 10, "vmin": -10}  # Default properties
    if plot_prop:
        plot_kwargs.update(plot_prop)

    # Plot the butterfly diagram
    mesh = ax.pcolormesh(ftime, y_vals, data, **plot_kwargs)

    # Set axis labels
    ax.set_xlabel("Time [Years]")
    ax.set_ylabel(y_label)

    # Add colorbar
    cbar_kwargs = {"label": "Magnetic Field Strength (G)"}
    if colorbar_prop:
        cbar_kwargs.update(colorbar_prop)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.01, **cbar_kwargs)

    # Add secondary axis for Carrington rotation number if provided
    if crn is not None:
        ax2 = ax.twiny()  # Create a twin x-axis
        ax2.set_xlim(ax.get_xlim())  # Match the x-limits with the primary axis

        # Map time to Carrington Rotation Number
        crn_ticks = np.linspace(ftime[0], ftime[-1], num=6)  # Define tick locations
        crn_labels = np.interp(crn_ticks, ftime, crn).astype(int)  # Interpolate CR numbers
        ax2.set_xticks(crn_ticks)
        ax2.set_xticklabels(crn_labels)
        ax2.set_xlabel("Carrington Rotation Number")

    return fig, ax



def plot_polar_field(ftime, north_pole, south_pole, shades=False, carrington_period=27.2753,
                     ax=None, plot_kwargs=None, label_suffix=""):
    """
    Plots the polar magnetic field strength over time, applying a moving average over one Carrington rotation.

    The function smooths the magnetic field strength at the north and south poles using a moving average
    with a window size equivalent to one Carrington rotation period. It then plots the smoothed values
    and optionally shades the regions below each curve for visualization.

    Parameters
    ----------
    label_suffix : str, optional
    ftime : ndarray
        1D array of time values (e.g., days or Carrington rotations).
    north_pole : ndarray
        1D array of magnetic field strength at the north pole.
    south_pole : ndarray
        1D array of magnetic field strength at the south pole.
    shades : bool, optional, default=False
        If True, fills the area below the plotted lines for better visualization.
    carrington_period : float, optional, default=27.2753
        Length of one Carrington rotation (in days) used for smoothing.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, a new figure is created.
    plot_kwargs : dict, optional
        Dictionary of keyword arguments passed to `ax.plot()` for customizing the plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The primary axes object of the plot.

    Example
    -------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> time = np.linspace(2000, 2025, 500)
    >>> north_pole = np.sin(2 * np.pi * time / 11) * 50 + np.random.randn(500) * 5
    >>> south_pole = -np.sin(2 * np.pi * time / 11) * 50 + np.random.randn(500) * 5
    >>> fig, ax = plot_polar_field(time, north_pole, south_pole, shades=True)
    >>> plt.show()
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure  # Get figure from the provided axes

    def moving_average(data, ftime, window_days):
        """Applies a moving average over a given time window in days."""
        window_size = np.searchsorted(ftime - ftime[0], window_days / 365.25)
        if window_size < 1:
            window_size = 1  # Ensure minimum valid window
        return np.convolve(data, np.ones(window_size) / window_size, mode="valid")

    # Apply moving average over one Carrington rotation
    north_pole_smooth = moving_average(north_pole, ftime, carrington_period)
    south_pole_smooth = moving_average(south_pole, ftime, carrington_period)

    # Adjust time array for valid moving average points
    time_smooth = ftime[len(ftime) - len(north_pole_smooth):]

    # Default plot properties
    plot_defaults = {"lw": 1.5, }
    if plot_kwargs:
        plot_defaults.update(plot_kwargs)

    # Plot smoothed north and south pole field strength
    ax.plot(time_smooth, north_pole_smooth, color="#B40426", label="North "+label_suffix, **plot_defaults)
    ax.plot(time_smooth, south_pole_smooth, color="#1F78B4", label="South "+label_suffix, **plot_defaults)

    # Fill area between lines for a smooth gradient effect
    if shades:
        ax.fill_between(time_smooth, south_pole_smooth, 0, color="#1F78B4", alpha=0.05)
        ax.fill_between(time_smooth, north_pole_smooth, 0, color="#B40426", alpha=0.05)

    # Add horizontal reference line at zero
    ax.axhline(0, color="gray", linestyle="dashed", linewidth=1)

    # Labels and Aesthetics
    ax.set_xlabel("Time [Years]")
    ax.set_ylabel("Polar Magnetic Field Strength (G)")
    ax.legend(frameon=True, handlelength=0.7)

    return fig, ax


class Visulalization:
    """
    A helper class to Visulaize the AFTmaps object.

    Parameters
    ----------
    aft : aft.AFTmaps
        AFTmaps object.

    Methods
    ---------
     plot_butterfly(self, ax=None, **kwargs):
        To Plot butterfly Diagram.

    """

    def __init__(self, aftdb, **kwargs):
        self.aftdb = aftdb
        extn = Path(aftdb).suffix[1:].upper()
        if extn=="JSON":
            self._data = pd.read_json(self.aftdb)
        elif extn=="ASDF":
            self._data = asdf.open("/Users/bjha/Data/AFT/AFTmapv2/SFTdb_20250311180144.asdf")
        else:
            NotImplementedError("Not Implemented Yet!!")
