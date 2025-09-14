import numpy as np
from bounded_minimize import bounded_minimize
from scipy.optimize import curve_fit
from scipy.stats import linregress
from uncertainties import ufloat
import uncertainties.umath as unp

def linear_function(x, a, b):
    return a * x + b

class Peak:
    def __init__(self, f_real, f_imag, peak_position, lower_continuum_edge, scaling=1):
        self.f_real = f_real
        self.f_imag = f_imag
        self.peak_position = peak_position
        self.lower_continuum_edge = lower_continuum_edge
        self.scaling = scaling
    
    def improved_peak_position(self, xtol=2e-12, offset = 0.2):
        offset_peak = offset * self.scaling
        search_bounds = (max(0, self.peak_position - offset_peak), min(self.peak_position + offset_peak, self.lower_continuum_edge - 0.01 * (self.lower_continuum_edge - self.peak_position)))
        
        result = bounded_minimize(self.f_imag, bounds=search_bounds, xtol=xtol)
        #self.peak_position = result["x"]
        return result
    
    def fit_real_part(self, range=0.001, begin_offset=1e-10, reversed=False, func=linear_function):
        """ If we fit
        ln(Re[G(z - z0)]) = a * ln(z - z0) + b
        
        and get a = -1, then the peak in the spectral function is a delta peak with weight e^b.
        This can be proven by means of the Kramers-Kronig relations. The factor 1/pi is absorbed in the definition of the spectral function.
        
        If the result is a=-2, then the peak is the derivative of a delta function.
        Such a peak does not have any weight, but a prefactor of e^b
        """
        begin_offset *= self.scaling
        range *= self.scaling
        lower_range = np.log(begin_offset)

        w_log = np.linspace(lower_range, np.log(begin_offset + range), 200)
        if reversed:
            w_usage = self.peak_position - np.exp(w_log)
        else:
            w_usage = self.peak_position + np.exp(w_log)
        
        # The absolute value is taken in case the real part is negative
        # This can be the case, depending on which side of the peak we are on
        # usually, if z>0 and if we are on the right side of the peak, real(data) > 0 and if we are left of the peaj real(data) < 0
        y_data = np.log(np.abs(self.f_real(w_usage)))
        return (linregress(w_log, y_data), w_log, y_data)

class PeakData:
    def __init__(self, position, weight, weight_error, slope):
        self.position = position
        self.weight = weight
        self.weight_error = weight_error
        self.slope = slope

    def __str__(self):
        return f"position = {self.position}\nweight = {self.weight}\nslope = {self.slope}"

def analyze_peak(f_real, f_imag, peak_position, lower_continuum_edge, 
                 peak_position_tol=1e-12, range=0.001, begin_offset=1e-10, 
                 reversed=False, improve_peak_position=True, peak_pos_range=0.01, plotter=None, scaling=1):
    """ Returns a PeakData object (peak position, peak weight, peak weight error)
    f_real                  Re[G(omega)]             We achieve better results if there is not imaginary shift here!
    f_imag                  Im[G(omega + i epsilon)] with some small positive constant epsilon
                                    Here, the imaginary shift is required to resolve delta peaks
    peak_position_tol       The tolerance for the peak position search
    range                   The range for the fit
    begin_offset            The offset for the range (for the phase peak at 0, this parameter should be larger than the default)
    reversed                Set True if the fit should be done on the left side of the peak
    improve_peak_position   By default the peak position is searched for numerical to improve precision
                            Set to False if this is undesired (e.g. you know the peak's position exactly)
                            If set to False, peak_position_tol does not do anything
    plotter                 If set, the fitted data and the fit are plotted using plotter.plot
    scaling                 Scales certain properties like the fit range, e.g., 1e-3 for conversion from meV to eV 
    """
    peak_position = peak_position.real # Assure that the peak position is a real number
    peak = Peak(f_real, f_imag, peak_position, lower_continuum_edge, scaling=scaling)
    peak_pos_value = np.copy(peak.peak_position)
    if improve_peak_position:
        peak_result = peak.improved_peak_position(xtol=peak_position_tol, offset=peak_pos_range)
        # only an issue if the difference is too large;
        if not peak_result["success"]:
            print("We might not have found the peak for data_folder!\nWe found ", peak_pos_value, " and\n", peak_result)
        # I am not quite sure, why this happens:
        # If the search range is too large, the algorithm might find a different mininum (that's logical)
        # However, if it find's such a "bad" minimum the return value seems to be always np.float(....)
        # If it finds the true mininum, the return value is np.array[...] of size 1 with the proper position
        # In both cases, the algorithm claims success....
        # So, since it seems to be good whenever we get an array as a result, we plainly use it.
        # The program dies if it only finds bad minima.
        if hasattr(peak_result["x"], "__len__"):
            omega_0 = peak_result["x"][0]
            peak.peak_position = omega_0
        else:
            peak_result = peak.improved_peak_position(xtol=peak_position_tol, offset=0.1 * peak_pos_range)
            omega_0 = peak_result["x"][0]
            peak.peak_position = omega_0
    else:
        omega_0 = peak_position
    
    result, w_log, y_data = peak.fit_real_part(range, begin_offset, reversed)
    if result.intercept > 36.8413614879:
        u_weight = ufloat(1e16, 1e16) # This is e^36.8413614879
    else:
        u_weight = unp.exp(ufloat(result.intercept, result.intercept_stderr))
    
    if plotter is not None:
        plotter.plot(w_log, y_data, label="Data")
        plotter.plot(w_log, linear_function(w_log, result.slope, result.intercept), label="Fit", linestyle="--")
        plotter.legend()
    return PeakData(omega_0, u_weight.nominal_value, u_weight.std_dev, ufloat(result.slope, result.stderr))