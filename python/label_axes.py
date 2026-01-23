import string
def label_axes(axes, x=0.05, y=0.95, va='top', ha='left', **kwargs):
    """Label axes in a grid of subplots.
    If axes is 2D, the default is 
    (a.1) (a.2) ...
    (b.1) (b.2) ...
    
    If axes is 1D, the default is
    (a) (b) (c) ...
    """
    
    if len(axes.shape) == 2:
        nrows, ncols = axes.shape
        for i in range(nrows):
            for j in range(ncols):
                label = f"({string.ascii_lowercase[i]}.{j+1})"
                axes[i, j].text(x, y, label, transform=axes[i, j].transAxes, va=va, ha=ha, **kwargs)
    elif len(axes.shape) == 1:
        ncols = axes.shape[0]
        for j in range(ncols):
            label = f"({string.ascii_lowercase[j]})"
            axes[j].text(x, y, label, transform=axes[j].transAxes, va=va, ha=ha, **kwargs)
    else:
        raise ValueError("Axes must be 1D or 2D array of matplotlib axes.")
    