import pandas as pd
import numpy as np
import continued_fraction_pandas as cf
from bounded_minimize import bounded_minimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.umath as unp

def linear_function(x, a, b):
    return a * x + b

class Peak:
    def __init__(self, f_real, f_imag, peak_position, lower_continuum_edge, imaginary_offset=1e-6):
        self.f_real = f_real
        self.f_imag = f_imag
        self.peak_position = peak_position
        self.lower_continuum_edge = lower_continuum_edge
        self.imaginary_offset = imaginary_offset
    
    def improved_peak_position(self, xtol=2e-12):
        offset_peak = 0.2
        search_bounds = (0 if self.peak_position - offset_peak < 0 else self.peak_position - offset_peak, 
                        self.lower_continuum_edge if self.peak_position + offset_peak > self.lower_continuum_edge else self.peak_position + offset_peak)
        
        result = bounded_minimize(self.f_imag, bounds=search_bounds, xtol=xtol)
        self.peak_position = result["x"]
        return result
    
    def fit_real_part(self, range=0.01, begin_offset=1e-10, reversed=False, func=linear_function):
        """ If we fit
        ln(Re[G(z - z0)]) = a * ln(z - z0) + b
        
        and get a = -1, then the peak in the spectral function is a delta peak with weight e^b.
        This can be proven by means of the Kramers-Kronig relations. The factor 1/pi is absorbed in the definition of the spectral function.
        
        If the result is a=-2, then the peak is the derivative of a delta function.
        Such a peak does not have any weight, but a prefactor of e^b
        """
        lower_range = np.log(begin_offset)

        w_log = np.linspace(lower_range, np.log(begin_offset + range), 2000, dtype=complex)
        w_log += (self.imaginary_offset * 1j)
        if reversed:
            w_usage = self.peak_position - np.exp(w_log)
        else:
            w_usage = self.peak_position + np.exp(w_log)
        
        # The absolute value is taken in case the real part is negative
        # This can be the case, depending on which side of the peak we are on
        # usually, if z>0 and if we are on the right side of the peak, real(data) > 0 and if we are left of the peaj real(data) < 0
        y_data = np.log(np.abs(self.f_real(w_usage)))
        self.popt, self.pcov = curve_fit(func, w_log, y_data)
        return self.popt, self.pcov, w_log, y_data

class PeakData:
    def __init__(self, position, weight, weight_error):
        self.position = position
        self.weight = weight
        self.weight_error = weight_error

def analyze_peak(f_real, f_imag, peak_position, lower_continuum_edge, imaginary_offset=1e-6, peak_position_tol=1e-12, range=0.001, begin_offset=1e-10, 
                 reversed=False, expected_slope=-1, plotter=None):
    """ Returns a PeakData object (peak position, peak weight, peak weight error)
    imaginary_offset is used to avoid the singularity at the peak position: z = omega + i * imaginary_offset
    peak_position_tol is the tolerance for the peak position search
    range is the range for the fit
    begin_offset is the offset for the range (for the phase peak at 0, this parameter should be larger than the default)
    reversed is True if the fit should be done on the left side of the peak
    expected_slope is the expected slope of the fit. If the slope is not as expected, a warning is printed.
    A slope of -1 indicates a delta peak with weight e^b
    """
    peak = Peak(f_real, f_imag, peak_position, lower_continuum_edge, imaginary_offset=imaginary_offset)
    peak_pos_value = np.copy(peak.peak_position)
    peak_result = peak.improved_peak_position(xtol=peak_position_tol)
    # only an issue if the difference is too large;
    if not peak_result["success"]:
        print("We might not have found the peak for data_folder!\nWe found ", peak_pos_value, " and\n", peak_result)

    popt, pcov, w_log, y_data = peak.fit_real_part(range, begin_offset, reversed)
    if abs(popt[0] - expected_slope) > 0.01: print("Warning!  Fit result: ", popt, "   But a slope of ", expected_slope, " was expected!")
    u_weight = unp.exp(ufloat(popt[1], np.sqrt(pcov[1][1])))
    
    if plotter is not None:
        plotter.plot(w_log, y_data, label="Data")
        plotter.plot(w_log, linear_function(w_log, *popt), label="Fit")
        plotter.legend()
    
    return PeakData(peak_result["x"], u_weight.nominal_value, u_weight.std_dev)