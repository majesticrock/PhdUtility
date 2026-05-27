import numpy as np
import string

def get(arr, i, j=0):
    arr = np.asarray(arr)
    if arr.ndim == 0:
        return arr.item()
    if arr.ndim == 1:
        if len(arr) == 0:
            return None
        return arr[i]
    if arr.ndim == 2:
        return arr[i, j]
    

def label_axes(axes, x=0.05, y=0.95, va='top', ha='left', special_labels=[], **kwargs):
    """Label axes in a grid of subplots.
    If axes is 2D, the default is 
    (a.1) (a.2) ...
    (b.1) (b.2) ...
    
    If axes is 1D, the default is
    (a) (b) (c) ...
    """
    axes = np.asarray(axes)
    if len(axes.shape) == 2:
        nrows, ncols = axes.shape
        for i in range(nrows):
            for j in range(ncols):
                extra = get(special_labels, i, j) if get(special_labels, i, j) is not None else ""
                label = f"({string.ascii_lowercase[i]}.{j+1}){extra}"
                axes[i, j].text(get(x, i, j), get(y, i, j), label, transform=axes[i, j].transAxes, va=va, ha=ha, **kwargs)
    elif len(axes.shape) == 1:
        ncols = axes.shape[0]
        for j in range(ncols):
            extra = get(special_labels, j) if get(special_labels, j) is not None else ""
            label = f"({string.ascii_lowercase[j]}){extra}"
            axes[j].text(get(x, j), get(y, j), label, transform=axes[j].transAxes, va=va, ha=ha, **kwargs)
    else:
        raise ValueError("Axes must be 1D or 2D array of matplotlib axes.")
    