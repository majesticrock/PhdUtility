import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NORM_FACTOR = -(1. / np.pi) 

class ContinuedFraction:
    def __init__(self, data_frame, z_squared=True, messages=True, ignore_first=5):
        self.z_squared = z_squared
        self.messages = messages
        self.ignore_first = ignore_first
        
        if z_squared:
            self.roots = np.array(data_frame.at["continuum_boundaries"])**2
        else:
            self.roots = np.array(data_frame.at["continuum_boundaries"])
        self.a_infinity = (self.roots[1] + self.roots[0]) * 0.5
        self.b_infinity = (self.roots[1] - self.roots[0]) * 0.25
        
        self.data = data_frame
        self.terminate_at = dict()

    def terminator(self, w_param):
        w = w_param**2
        p = w - self.a_infinity
        q = 4 * self.b_infinity**2
        root = np.sqrt(np.real(p**2 - q), dtype=complex)
        
        if hasattr(w, '__len__'):
            return_arr = np.zeros(len(w), dtype=complex)
            for i in range(0, len(w)):
                if w_param[i].real > 0:
                    if w[i].real > self.roots[0]:
                        return_arr[i] = (p[i] - root[i]) / (2. * self.b_infinity**2)
                    else:
                        return_arr[i] = (p[i] + root[i]) / (2. * self.b_infinity**2)
                else:
                    if w[i].real > self.roots[1]:
                        return_arr[i] = (p[i] - root[i]) / (2. * self.b_infinity**2)
                    else:
                        return_arr[i] = (p[i] + root[i]) / (2. * self.b_infinity**2)
        else:
            if w_param.real > 0:
                if w.real > self.roots[0]:
                    return (p - root) / (2. * self.b_infinity**2)
                else:
                    return (p + root) / (2. * self.b_infinity**2)
            else:
                if w.real > self.roots[1]:
                    return (p - root) / (2. * self.b_infinity**2)
                else:
                    return (p + root) / (2. * self.b_infinity**2)
        return return_arr
    
    def __coeffs_A__(self, name, index):
        return self.data[f"resolvents.{name}"][index]["a_i"]
    
    def __coeffs_B__(self, name, index):
        return self.data[f"resolvents.{name}"][index]["b_i"]
    
    def find_termination_depth(self, name, index):
        A = self.__coeffs_A__(name, index)
        B = self.__coeffs_B__(name, index)
        
        deviation_from_infinity = np.zeros(len(A) - 1)
        for i in range(0, len(A) - 1):
            deviation_from_infinity[i] = abs((A[i] - self.a_infinity) / self.a_infinity) + abs((np.sqrt(B[i + 1]) - self.b_infinity) / self.b_infinity)

        best_approx = np.argmin(deviation_from_infinity[self.ignore_first:]) + self.ignore_first
        #best_approx = 45
        self.terminate_at[name][index] = len(A) - best_approx
        if self.messages: 
            print("Terminating at i =", best_approx)
    
    def __termination_depth_if_required__(self, name, index):
        if not name in self.terminate_at:
            self.terminate_at[name] = np.zeros(len(self.data[f"resolvents.{name}"]), dtype=int)
        if self.terminate_at[name][index] == 0:
            self.find_termination_depth(name, index)
    
    def continued_fraction(self, w_param, name, index=0, withTerminator=True):
        self.__termination_depth_if_required__(name, index)
        A = self.__coeffs_A__(name, index)
        B = self.__coeffs_B__(name, index)
        w = w_param**2
        if withTerminator:
            G = w - A[len(A) - self.terminate_at[name][index]] - B[len(B) - self.terminate_at[name][index]] * self.terminator(w_param.real)
        else:
            G = w - A[len(A) - self.terminate_at[name][index]]
        for j in range(len(A) - self.terminate_at[name][index] - 1, -1, -1):
            G = w - A[j] - B[j + 1] / G
        return B[0] / G
    
    def spectral_density(self, w_param, name, index=0, withTerminator = True):
        return NORM_FACTOR * self.continued_fraction(w_param, name, index, withTerminator).imag
    
    def denominator(self, w_param, name, index=0, withTerminator=True):
        self.__termination_depth_if_required__(name, index)
        A = self.__coeffs_A__(name, index)
        B = self.__coeffs_B__(name, index)
        w = w_param**2
        if withTerminator:
            G = w - A[len(A) - self.terminate_at[name][index]] - B[len(B) - self.terminate_at[name][index]] * self.terminator(w_param.real)
        else:
            G = w - A[len(A) - self.terminate_at[name][index]]
        for j in range(len(A) - self.terminate_at[name][index] - 1, -1, -1):
            G = w - A[j] - B[j + 1] / G
        return G / B[0]
    
    def mark_continuum(self, axes=None, scale_factor=1., label="Continuum"):
        if label is not None:
            args = {"alpha" : 0.333, "color": "grey", "label" : label}
        else:
            args = {"alpha" : 0.333, "color": "grey"}
            
        if axes is None:
            plotter = plt.axvspan
        else:
            plotter = axes.axvspan
            
        if self.z_squared:
            plotter(scale_factor * np.sqrt(self.roots[0]), scale_factor * np.sqrt(self.roots[1]), **args)
        else:
            plotter(scale_factor * self.roots[0], scale_factor * self.roots[1], **args)
    
    def continuum_edges(self):
        if not self.z_squared:
            return self.roots
        else:
            return np.sqrt(self.roots)