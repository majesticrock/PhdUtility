import numpy as np

class SingleParticle:
    def __init__(self, gap_parameters):
        self.gap_squared = gap_parameters[0]**2 + gap_parameters[1]**2 + gap_parameters[2]**2
        self.cos_factor = -(2. + 0.5 * (gap_parameters[6] + gap_parameters[7]))
    
    def bare_dispersion(self, _kx, _ky):
        X, Y = np.meshgrid(_kx, _ky)
        return self.cos_factor * (np.cos(np.pi * X) + np.cos(np.pi * Y))
    
    def dispersion(self, _kx, _ky):
        return np.sqrt(self.gap_squared + self.bare_dispersion(_kx, _ky)**2)

def continuum_bounds(_dispersion, _kx, _ky, N_k):
    qx = np.linspace(0, 1, N_k)
    qy = np.linspace(0, 1, N_k)
    
    diff_omega = _dispersion(qx, qy) + _dispersion(_kx - qx, _ky - qy)
    return np.min(diff_omega), np.max(diff_omega)