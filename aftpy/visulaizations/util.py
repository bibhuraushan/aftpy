import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean

cmap_balance = cmocean.cm.balance

def cmap_YlGrBl():
    """
    Creates a custom colormap where negative values are blue,
    positive values are yellow, and zero is light gray.

    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
        The generated colormap.
    """
    colors = [
        (0, "blue"),  # Blue at negative max
        (0.5, "#848884"),  # Light gray at zero
        (1, "yellow")  # Yellow at positive max
    ]
    cmap = mcolors.LinearSegmentedColormap.from_list("butterfly_cmap", colors)

    # cmap = mcolors.LinearSegmentedColormap.from_list("blue_gray_yellow", colors)
    return cmap

