import numpy as np
import pandas as pd

NORM_FACTOR = -(1. / np.pi) 

class ContinuedFraction:
    def __init__(self, data_frame, z_squared=True, messages=True, ingore_first=40):
        self.z_squared = z_squared
        self.messages = messages
        
        if z_squared:
            self.roots = np.array(data_frame.at["Continuum Boundaries"])**2
        else:
            self.roots = np.array(data_frame.at["Continuum Boundaries"])
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
    
    def continued_fraction(self, w_param, name, withTerminator=True):
        w = w_param**2
        if withTerminator:
            G = w - self.A[len(self.A) - self.terminate_at] - self.B[len(self.B) - self.terminate_at] * self.terminator(w_param.real)
        else:
            G = w - self.A[len(self.A) - self.terminate_at]
        for j in range(len(self.A) - self.terminate_at - 1, -1, -1):
            G = w - self.A[j] - self.B[j + 1] / G

        return self.B[0] / G