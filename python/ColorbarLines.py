import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from collections.abc import Iterable
import matplotlib as mpl
import numpy as np

class ColorbarLines:
    def __init__(self, fig, cbar_axes, n, vmin=0, vmax=1, vcenter=None, norm=None, cmap='viridis', cbar_label=""):
        """
        fig:        figure to plot on
        cbar_axes:  axes the colorbar will steal figure space from
        n:          number of colorbar levels
        vmin/vmax:  mininum and maximum colorbar level
        vcenter:    middle colorbar level. If set, the default will use a TwoSlopeNorm
        norm:       a custom norm, if desired. If set, vmin, vmax, and vcenter loose their meaning
        cmap:       colormap to draw the colors from
        cbar_label: label for the colorbar
        
        The colorbar itself can be accessed as a member named self.cbar
        """
        self.cmap = mpl.colormaps[cmap] if isinstance(cmap, str) else cmap
        
        self.colors = self.cmap(np.linspace(0, 1, n))
        listed_cmap = ListedColormap(self.colors)
        if norm is None:
            if vcenter is None:
                norm = plt.Normalize(vmin=vmin, vmax=vmax)
            else:
                norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=listed_cmap, norm=norm)
        self.cbar_space = norm.inverse(np.linspace(0, 1, n))
        
        if cbar_axes is not None:
            self.cbar = fig.colorbar(sm, ax=cbar_axes, label=cbar_label)
        
    def plot(self, xs, ys, color_values, ax, **kwargs):
        """
        Plots (x_i, y_i) and determines the colors based on color_values.
        The idea is, that the arrays x_i and y_i correspond to the datapoint in color_values[i].
        
        xs:             x-values - either a 2D array of the same shape as ys or a 1D array of the same shapes as ys[i]
        ys:             y-values
        color_values:   values that determine the lines' positions on the colorscale
        ax:             axis to plot on
        kwargs:         kwargs forwarded to the plot function
        """
        if isinstance(xs[0], Iterable):
            assert len(xs) == len(ys)
        assert len(color_values) == len(ys)
        length = len(ys)  
        
        for i in range(length):
            idx = (np.abs(self.cbar_space - color_values[i])).argmin()
            if isinstance(xs[0], Iterable):
                ax.plot(xs[i], ys[i], c=self.colors[idx], **kwargs)
            else:
                ax.plot(xs, ys[i], c=self.colors[idx], **kwargs)


if __name__ == '__main__':
    fig, axes = plt.subplots(ncols=3, figsize=(12.4, 4.8))
    x = np.linspace(-3, 3)
    A = np.linspace(-10, 20, 50)
    
    cbar_plotter = ColorbarLines(fig, axes, vmin=-10, vmax=20, n=500, cmap='seismic', cbar_label="a", vcenter=0)
    
    def f(x, a):
        return a * np.cos(x)
    
    ys = np.array([ f(x, a) for a in A ])
    cbar_plotter.plot(x, ys, A, axes[0])
    
    T = np.linspace(-np.pi, np.pi, len(A))
    xs = np.array([ x + t for t in T ])
    ys = np.array([ f(__x, a) for __x, a in zip(xs, A) ])
    cbar_plotter.plot(xs, ys, A, axes[1])
    
    A2 = np.linspace(-3, 10, 25)
    ys = np.array([ f(x, a) for a in A2 ])
    cbar_plotter.plot(x, ys, A2, axes[2])
    
    plt.show()