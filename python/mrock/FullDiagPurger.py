import numpy as np

PEAK_DIFF_TOL = 0.004

def __fill_temps__(eigenvalues, weights, first_eigenvectors, continuum_edge):
    temp_evs = []
    temp_weights = []
    temp_vectors = []
    for ev, weight, vector in zip(eigenvalues, weights, first_eigenvectors):
        # if the weight is less than machine epsilon, it's an artifact
        if ev >= continuum_edge or weight < 1e-16:
            continue
        if len(temp_evs) == 0:
            temp_evs.append(ev)
            temp_weights.append(weight)
            temp_vectors.append(vector)
        else:
            if weight < 1e-8:
                continue
            if abs(ev - temp_evs[-1]) > PEAK_DIFF_TOL:
                temp_evs.append(ev)
                temp_weights.append(weight)
                temp_vectors.append(vector)
            elif weight > temp_weights[-1]:
                temp_evs[-1] = ev
                temp_weights[-1] = weight
                temp_vectors[-1] = vector
    return temp_evs, temp_weights, temp_vectors

class FullDiagPurger:
    def __init__(self, system_data, epsilon_space):
        temp_amplitude_evs, temp_amplitude_weights, temp_amplitude_vectors = __fill_temps__(
            system_data["amplitude.eigenvalues"], 
            system_data["amplitude.weights"][0], 
            system_data["amplitude.first_eigenvectors"],
            system_data["continuum_boundaries"][0]
        )
        temp_phase_evs, temp_phase_weights, temp_phase_vectors = __fill_temps__(
            system_data["phase.eigenvalues"], 
            system_data["phase.weights"][0], 
            system_data["phase.first_eigenvectors"],
            system_data["continuum_boundaries"][0]
        )
        self.epsilon_space = epsilon_space
        self.N = len(epsilon_space)
        
        doublets = []
        for i, a_ev in enumerate(temp_amplitude_evs):
            for j, p_ev in enumerate(temp_phase_evs):
                if abs(a_ev - p_ev) < PEAK_DIFF_TOL:
                    doublets.append((i, j))

        correct_amplitude_indices = [ i for i in range(len(temp_amplitude_evs)) if all(i != dt[0] for dt in doublets) ]
        correct_phase_indices = [ j for j in range(len(temp_phase_evs)) if all(j != dt[1] for dt in doublets) ]
        glimmering_amplitude_indices = []
        glimmering_phase_indices = []
        
        for doublet in doublets:
            a_i, p_i = doublet
            if temp_amplitude_weights[a_i] < temp_phase_weights[p_i]:
                correct_phase_indices.append(p_i)
                glimmering_amplitude_indices.append(a_i)
            else:
                correct_amplitude_indices.append(a_i)
                glimmering_phase_indices.append(p_i)
        
        correct_amplitude_indices.sort()
        correct_phase_indices.sort()
        
        self.amplitude_eigenvalues  = np.array([ temp_amplitude_evs[i]     for i in correct_amplitude_indices ])
        self.amplitude_weights      = np.array([ temp_amplitude_weights[i] for i in correct_amplitude_indices ])
        self.amplitude_eigenvectors = np.array([ temp_amplitude_vectors[i] for i in correct_amplitude_indices ])
        
        self.phase_eigenvalues  = np.array([ temp_phase_evs[i]     for i in correct_phase_indices ])
        self.phase_weights      = np.array([ temp_phase_weights[i] for i in correct_phase_indices ])
        self.phase_eigenvectors = np.array([ temp_phase_vectors[i] for i in correct_phase_indices ])
        
        for i in range(len(self.amplitude_eigenvectors)):
            if np.sum(self.amplitude_eigenvectors[i][:self.N]) < 0:
                self.amplitude_eigenvectors[i] *= -1
        for i in range(len(self.phase_eigenvectors)):
            if np.sum(self.phase_eigenvectors[i]) < 0:
                self.phase_eigenvectors[i] *= -1
        
        
                
        glimmering_amplitude_indices.sort()
        glimmering_phase_indices.sort()
        
        self.glimmering_amplitude_eigenvalues  = np.array([ temp_amplitude_evs[i]     for i in glimmering_amplitude_indices ])
        self.glimmering_amplitude_weights      = np.array([ temp_amplitude_weights[i] for i in glimmering_amplitude_indices ])
        self.glimmering_amplitude_eigenvectors = np.array([ temp_amplitude_vectors[i] for i in glimmering_amplitude_indices ])
        
        self.glimmering_phase_eigenvalues  = np.array([ temp_phase_evs[i]     for i in glimmering_phase_indices ])
        self.glimmering_phase_weights      = np.array([ temp_phase_weights[i] for i in glimmering_phase_indices ]) 
        self.glimmering_phase_eigenvectors = np.array([ temp_phase_vectors[i] for i in glimmering_phase_indices ])
        
        for i in range(len(self.glimmering_amplitude_eigenvectors)):
            if np.sum(self.glimmering_amplitude_eigenvectors[i][:self.N]) < 0:
                self.glimmering_amplitude_eigenvectors[i] *= -1
        for i in range(len(self.glimmering_phase_eigenvectors)):
            if np.sum(self.glimmering_phase_eigenvectors[i]) < 0:
                self.glimmering_phase_eigenvectors[i] *= -1
    
    def plot_line(self, ax, y, combined_norm=True, **plot_kwargs):
        if len(y) == 2 * self.N:
            if combined_norm:
                norm = 1. / np.max(np.abs(y))
                ax[0].plot(self.epsilon_space, y[:self.N] * norm, **plot_kwargs)
                ax[1].plot(self.epsilon_space, y[self.N:] * norm, **plot_kwargs)
                return norm
            else:
                norm1 = 1. / np.max(np.abs(y[:self.N]))
                norm2 = 1. / np.max(np.abs(y[self.N:]))
                ax[0].plot(self.epsilon_space, y[:self.N] * norm1, **plot_kwargs)
                ax[1].plot(self.epsilon_space, y[self.N:] * norm2, **plot_kwargs)
                return norm1 / norm2
        else:
            norm = 1. / np.max(np.abs(y))
            ax.plot(self.epsilon_space, y * norm, **plot_kwargs)
            return norm
    
    ##############################################################
    def __plot_impl__(self, ax, eigenvectors, eigenvalues, which, label_energy, combined_norm=True, **plot_kwargs):
        if which == 'all':
            for i in range(len(eigenvectors)):
                if label_energy:
                    plot_kwargs["label"] = f"{i+1} ({eigenvalues[i]:.3f})"
                else:
                    plot_kwargs["label"] = f"{i+1}"
                norm = self.plot_line(ax, eigenvectors[i], combined_norm, **plot_kwargs)
                if combined_norm == False:
                    print(f"Multiplied A_eps for i={i} by a factor of {norm}")
        else:
            if len(eigenvectors) <= which:
                print(f"Warning: Tried to plot eigenvector {which}, but only {len(eigenvectors)} available.")
            else:
                if label_energy:
                    plot_kwargs['label'] = f"{eigenvalues[which]:.3f}"
                norm = self.plot_line(ax, eigenvectors[which], combined_norm, **plot_kwargs)
                if combined_norm == False:
                    print(f"Multiplied A_eps for i={which} by a factor of {norm}")
    
    ##############################################################
    def plot_phase(self, ax, which='all', label_energy=False, **plot_kwargs):
        self.__plot_impl__(ax, self.phase_eigenvectors, self.phase_eigenvalues, which, label_energy, **plot_kwargs)
    
    ##############################################################
    def plot_glimmering_phase(self, ax, which='all', label_energy=False, **plot_kwargs):
        self.__plot_impl__(ax, self.glimmering_phase_eigenvectors, self.glimmering_phase_eigenvalues, which, label_energy, **plot_kwargs)
    
    ##############################################################
    def plot_amplitude(self, axes, which='all', label_energy=False, combined_norm=True, **plot_kwargs):
        self.__plot_impl__(axes, self.amplitude_eigenvectors, self.amplitude_eigenvalues, which, 
                           label_energy, combined_norm, **plot_kwargs)
    
    ##############################################################
    def plot_glimmering_amplitude(self, axes, which='all', label_energy=False, combined_norm=True, **plot_kwargs):
        self.__plot_impl__(axes, self.glimmering_amplitude_eigenvectors, self.glimmering_amplitude_eigenvalues, which, 
                           label_energy, combined_norm, **plot_kwargs)
        
    # returns a tuple (sum_k alpha_k^2, sum_k nu_k^2)
    def integral_amplitude(self, which):
        if which < len(self.amplitude_eigenvectors):
            pc = np.sum(self.amplitude_eigenvectors[which][:self.N]**2)
            n  = np.sum(self.amplitude_eigenvectors[which][self.N:]**2)
            #print(f"Pair creation integral: {pc}\nOccupation integral: {n}")
            #print(f"Combined: {pc + n}")
            return pc, n
        else:
            raise ValueError(f"Requested integral for amplitude mode {which}, but only {len(self.amplitude_eigenvectors)} available.")
        
    def to_dict(self, suffix):
        return { 
                f"energies_higgs_{suffix}" : self.amplitude_eigenvalues,
                f"energies_phase_{suffix}" : self.phase_eigenvalues,
                f"amplitudes_higgs_{suffix}" : self.amplitude_eigenvectors,
                f"amplitudes_phase_{suffix}" : self.phase_eigenvectors
            }