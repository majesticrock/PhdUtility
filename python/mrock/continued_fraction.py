import numpy as np
import matplotlib.pyplot as plt

from . import cpp_continued_fraction as ccf

NORM_FACTOR = -(1. / np.pi) 

class ContinuedFraction:
    def __init__(self, data_frame, messages=True, ignore_first=5, ignore_last=80):
        self.messages = messages
        self.ignore_first = ignore_first
        self.ignore_last = ignore_last
        
        boundaries = np.asarray(data_frame.at["continuum_boundaries"], dtype=float)
        if boundaries.shape != (2,):
            raise ValueError("continuum_boundaries must contain exactly two values.")
        if boundaries[0] >= boundaries[1]:
            raise ValueError("continuum_boundaries must be ordered as [lower, upper].")
        self.continuum_boundaries_squared = boundaries**2

        self.a_infinity = (self.continuum_boundaries_squared[1] + self.continuum_boundaries_squared[0]) * 0.5
        self.b_infinity = (self.continuum_boundaries_squared[1] - self.continuum_boundaries_squared[0]) * 0.25
        
        self.data = data_frame
        self.termination_resolvent_index = dict()

    def terminator(self, w_param, resolvent_name, resolvent_index=0, with_terminator=True):
        data = self._get_cf_data(resolvent_name, resolvent_index, with_terminator)
        w_param = np.asarray(w_param, dtype=np.float64)
        return ccf.terminator(w_param, data)
    
    def _coeffs_A(self, resolvent_name, resolvent_index):
        return np.asarray(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["a_i"])
    
    def _coeffs_B(self, resolvent_name, resolvent_index):
        return np.asarray(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["b_i"])
    
    def _len_A(self, resolvent_name, resolvent_index):
        return len(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["a_i"])
    
    def find_termination_depth(self, resolvent_name, resolvent_index):
        A = np.asarray(self._coeffs_A(resolvent_name, resolvent_index), dtype=np.float64)
        B = np.asarray(self._coeffs_B(resolvent_name, resolvent_index), dtype=np.float64)
        
        ignore_last = min(self.ignore_last, len(A) - 1)            
        deviation_from_infinity = np.zeros(self.ignore_last)
        for i in range(0, ignore_last):
            deviation_from_infinity[i] = abs((A[i] - self.a_infinity) / self.a_infinity) + abs((np.sqrt(B[i + 1]) - self.b_infinity) / self.b_infinity)

        if self.ignore_first >= ignore_last:
            raise ValueError(
                f"ignore_first={self.ignore_first} must be smaller than ignore_last={ignore_last}."
            )
        best_approx = np.argmin(deviation_from_infinity[self.ignore_first:]) + self.ignore_first

        self.termination_resolvent_index[resolvent_name][resolvent_index] = best_approx
        if self.messages: 
            print("Terminating at i =", len(A) - best_approx)
    
    def _termination_depth_if_required(self, resolvent_name, resolvent_index):
        if not resolvent_name in self.termination_resolvent_index:
            self.termination_resolvent_index[resolvent_name] = np.zeros(len(self.data[f"resolvents.{resolvent_name}"]), dtype=int)
        if self.termination_resolvent_index[resolvent_name][resolvent_index] == 0:
            self.find_termination_depth(resolvent_name, resolvent_index)
    
    def _get_cf_data(self, resolvent_name, resolvent_index, with_terminator):
        self._termination_depth_if_required(resolvent_name, resolvent_index)
        return ccf.ContinuedFractionData(self.a_infinity, self.b_infinity**2, np.asarray(self.continuum_boundaries_squared), 
                                         np.asarray(self._coeffs_A(resolvent_name, resolvent_index)), np.asarray(self._coeffs_B(resolvent_name, resolvent_index)), 
                                         self.termination_resolvent_index[resolvent_name][resolvent_index], with_terminator)
    
    def continued_fraction(self, w_param, resolvent_name, resolvent_index=0, with_terminator=True):
        scalar_input = np.isscalar(w_param)
        w_param = np.atleast_1d(np.asarray(w_param, dtype=np.complex128))
        result = ccf.continued_fraction(w_param, self._get_cf_data(resolvent_name, resolvent_index, with_terminator))
        return result[0] if scalar_input else result
    
    def continued_fraction_varied_depth(self, w_param, resolvent_name, shift_range, resolvent_index=0, with_terminator=True):
        scalar_input = np.isscalar(w_param)
        w_param = np.atleast_1d(np.asarray(w_param, dtype=np.complex128))
        result = ccf.continued_fraction_varied_depth(w_param, self._get_cf_data(resolvent_name, resolvent_index, with_terminator), shift_range)
        return result[0] if scalar_input else result

    def spectral_density(self, w_param, resolvent_name, resolvent_index=0, with_terminator=True):
        return NORM_FACTOR * self.continued_fraction(w_param, resolvent_name, resolvent_index, with_terminator).imag
    
    def spectral_density_varied_depth(self, w_param, resolvent_name, shift_range, resolvent_index=0, with_terminator=True):
        return NORM_FACTOR * self.continued_fraction_varied_depth(w_param, resolvent_name, shift_range, resolvent_index, with_terminator).imag
    
    def real_part(self, w_param, resolvent_name, resolvent_index=0, with_terminator=True):
        return self.continued_fraction(w_param, resolvent_name, resolvent_index, with_terminator).real
    
    def denominator(self, w_param, resolvent_name, resolvent_index=0, with_terminator=True):
        scalar_input = np.isscalar(w_param)
        w_param = np.atleast_1d(np.asarray(w_param, dtype=np.complex128))
        result = ccf.denominator(w_param, self._get_cf_data(resolvent_name, resolvent_index, with_terminator))
        return result[0] if scalar_input else result
    
    # Does not really work for Goldstone peaks
    def classify_bound_states(self, resolvent_name, resolvent_index=0, n_scan=1000, weight_domega=1e-8, tolerance_bits=48, max_iter=200, is_phase_peak=None):
        result = ccf.classify_bound_states(self._get_cf_data(resolvent_name, resolvent_index, with_terminator=True), n_scan, weight_domega, tolerance_bits, max_iter)
        if is_phase_peak is not None:
            result = [ item for item in result if not is_phase_peak(item[0]) ]
        return result
                
    
    def mark_continuum(self, ax=None, scale_factor=1., label="Continuum"):
        if label is not None:
            args = {"alpha" : 0.333, "color": "grey", "label" : label}
        else:
            args = {"alpha" : 0.333, "color": "grey"}
            
        if ax is None:
            plotter = plt.axvspan
        else:
            plotter = ax.axvspan
            
        plotter(scale_factor * np.sqrt(self.continuum_boundaries_squared[0]), scale_factor * np.sqrt(self.continuum_boundaries_squared[1]), **args)
    
    def continuum_edges(self):
        return np.sqrt(self.continuum_boundaries_squared)