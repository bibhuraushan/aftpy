import numpy as np
import pylab as p
from matplotlib import pyplot as plt


class Visulalization:
    def __init__(self, aftobj, **kwargs):
        self.aftobj = aftobj

    def plot_butterfly(self, ax=None, **kwargs):
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
