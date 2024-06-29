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
