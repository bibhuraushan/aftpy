import numpy as np
import pylab as p
from matplotlib import pyplot as plt
from matplotlib.axis import Axis


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

    def __init__(self, aftobj, **kwargs):
        self.aftobj = aftobj

    def plot_butterfly(self, ax: Axis = None, **kwargs) -> Axis:
        """
        Function To plot butterfly Diagram for giveb AFTmaps object.

        Parameters
        ----------
        ax: matplotlib Axes
            Axis to plot butterfly Diagram for.
        kwargs: dict
            Keyword Args to pass to matplotlib Axes.
        """
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            ax.set(xlabel="Time [Years]", ylabel="Latitude [degree]")
            bimage = self.aftobj.cravgmap()
            bimage = np.clip(bimage, -10, 10)
            _extent = (self.aftobj.crn.min(), self.aftobj.crn.max() + 1, -90, 90)
            ax.imshow(bimage, origin="lower", cmap="bwr",
                      aspect="auto", vmax=10, vmin=-10,
                      extent=_extent)
            p.show()
        else:
            raise NotImplementedError(f"Butterfly diagram not implemented for {ax}")

    def plot_polarfiled(self, ax: Axis = None, **kwargs):
        """

        Parameters
        ----------
        ax
        kwargs
        """
        raise NotImplementedError(f"Polar diagram not implemented for {ax}")


def plot_butterfly(data, time, crn=None, sinlat=True,
                   ax=None, plot_prop=None, colorbar_prop=None):
    """
    Plots the magnetic butterfly diagram with an optional secondary axis for Carrington rotation numbers.

    Parameters
    ----------
    data : ndarray
        2D array representing the magnetic field data (latitude vs. time).
    time : ndarray
        1D array of time values corresponding to the x-axis.
    crn : ndarray, optional
        1D array of Carrington rotation numbers (same length as `time`).
    sinlat : bool, optional
        If True, plot against sine(latitude); otherwise, use latitude directly.
    ax : matplotlib.axes._subplots.Axes, optional
        Axes object to plot on. If None, creates a new figure.
    plot_prop : dict, optional
        Dictionary of properties passed to `ax.pcolormesh()`.
    colorbar_prop : dict, optional
        Dictionary of properties passed to `fig.colorbar()`.

    Returns
    -------
    fig, ax
        Figure and Axes objects of the plot.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure  # Get figure from the provided axes

    cmap = create_blue_gray_yellow_colormap()

    # Define latitude range
    latitudes = np.linspace(-90, 90, data.shape[0])
    if sinlat:
        y_vals = np.sin(np.radians(latitudes))
        y_label = "Sine Latitude"
    else:
        y_vals = latitudes
        y_label = "Latitude (°)"

    # Default plotting properties
    plot_kwargs = {"cmap": butterfly_cmap, "vmax":10, "vmin":-10}  # Default properties
    if plot_prop:
        plot_kwargs.update(plot_prop)

    # Plot the butterfly diagram
    mesh = ax.pcolormesh(time, y_vals, data, **plot_kwargs)
    # Set axis labels
    ax.set_xlabel("Time [Years]")
    ax.set_ylabel(y_label)
    # ax.set_title("Magnetic Butterfly Diagram")

    # Add colorbar
    cbar_kwargs = {"label": "Magnetic Field Strength (G)"}
    if colorbar_prop:
        cbar_kwargs.update(colorbar_prop)
    cbar = fig.colorbar(mesh, ax=ax, **cbar_kwargs)

    # Add secondary axis for Carrington rotation number if provided
    if crn is not None:
        ax2 = ax.twiny()  # Create a twin x-axis
        ax2.set_xlim(ax.get_xlim())  # Match the x-limits with the primary axis

        # Map time to Carrington Rotation Number
        crn_ticks = np.linspace(time[0], time[-1], num=6)  # Define tick locations
        crn_labels = np.interp(crn_ticks, time, crn).astype(int)  # Interpolate CR numbers
        ax2.set_xticks(crn_ticks)
        ax2.set_xticklabels(crn_labels)
        ax2.set_xlabel("Carrington Rotation Number")

    return fig, ax