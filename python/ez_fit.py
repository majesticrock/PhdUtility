import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Creates a numpy linspace with the bounds of coarse_space
def ez_lin_space(coarse_space):
    return np.linspace(min(coarse_space), max(coarse_space), 200)

# Linearly fits to a data set and plots the fit, returns the best fit parameters
def ez_linear_fit(x_data, y_data, plotter=plt, x_space=None, **plotter_kwargs):
    def func(x, a, b):
        return a*x + b
    popt, pcov = curve_fit(func, x_data, y_data)
    
    if not hasattr(x_space, "__len__"):
        x_space = ez_lin_space(x_data)
    plotter.plot(x_space, func(x_space, *popt), **plotter_kwargs)
    
    return popt, pcov

def ez_general_fit(x_data, y_data, func, plotter=plt, x_space=None, **plotter_kwargs):
    popt, pcov = curve_fit(func, x_data, y_data)
    
    if not hasattr(x_space, "__len__"):
        x_space = ez_lin_space(x_data)
    plotter.plot(x_space, func(x_space, *popt), **plotter_kwargs)
    
    return popt, pcov