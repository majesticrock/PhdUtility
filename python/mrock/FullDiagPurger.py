"""Utilities for filtering and plotting the operator amplitudes computed via
the iEoM library.

This module provides tools for post-processing eigenvalues, spectral weights,
and eigenvectors obtained from full diagonalization data. It removes likely
artifacts, separates amplitude and phase modes, identifies near-degenerate
``doublets``, and stores weaker companion peaks as ``glimmering`` modes.

The main entry point is :class:`FullDiagPurger`.
"""

import numpy as np

PEAK_DIFF_TOL = 0.004


def _fill_temporaries(eigenvalues, weights, first_eigenvectors, continuum_edge):
    """Filter and merge candidate bound-state peaks below the continuum edge.

    This helper iterates through eigenvalues, spectral weights, and first
    eigenvectors and keeps only physically relevant candidate peaks. Peaks above
    the continuum edge or with negligible weight are discarded. Nearby peaks
    closer than :data:`PEAK_DIFF_TOL` are treated as the same peak, and only the
    one with the larger weight is retained.

    Parameters
    ----------
    eigenvalues : array-like
        Eigenvalues to inspect.
    weights : array-like
        Spectral weights associated with ``eigenvalues``.
    first_eigenvectors : array-like
        First eigenvectors associated with ``eigenvalues``.
    continuum_edge : float
        Lower continuum boundary. Eigenvalues greater than or equal to this
        value are discarded.

    Returns
    -------
    tuple[list, list, list]
        Tuple ``(temp_evs, temp_weights, temp_vectors)`` containing filtered
        eigenvalues, weights, and eigenvectors.

    Notes
    -----
    Weights below ``1e-16`` are considered numerical artifacts and discarded.
    After at least one peak has been accepted, additional peaks with weights
    below ``1e-8`` are skipped.
    """
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
    """Post-process the operator amplitudes of the amplitude and phase modes.

    The class filters raw full-diagonalization output, separates amplitude and
    phase modes, identifies nearly degenerate amplitude-phase pairs, and stores
    weaker members of such pairs as ``glimmering`` modes.
    
    These ``glimmering`` modes appear because any energy eigenvalue is present in
    both channels. If the system exhibits particle-hole symmetry, the amplitude modes'
    weights in the phase channel vanishes, and vice versa. If the system does not
    exhibit this symmetry, the weights in the ``wrong`` channel are nonzero, but 
    typically smaller than the weights in the ``correct`` channel.
    This fact is used to assign any given mode to either the amplitude or phase 
    channel, and to store the weaker companion mode.

    Parameters
    ----------
    system_data : mapping
        Dictionary-like object containing full-diagonalization data. Expected
        keys include:

        - ``"amplitude.eigenvalues"``
        - ``"amplitude.weights"``
        - ``"amplitude.first_eigenvectors"``
        - ``"phase.eigenvalues"``
        - ``"phase.weights"``
        - ``"phase.first_eigenvectors"``
        - ``"continuum_boundaries"``

    epsilon_space : array-like
        Energy or momentum grid used for plotting eigenvector components.

    Attributes
    ----------
    epsilon_space : array-like
        Grid associated with the eigenvector components.
    N : int
        Length of ``epsilon_space``.
    amplitude_eigenvalues : numpy.ndarray
        Filtered amplitude-mode eigenvalues.
    amplitude_weights : numpy.ndarray
        Filtered amplitude-mode weights.
    amplitude_eigenvectors : numpy.ndarray
        Filtered amplitude-mode eigenvectors.
    phase_eigenvalues : numpy.ndarray
        Filtered phase-mode eigenvalues.
    phase_weights : numpy.ndarray
        Filtered phase-mode weights.
    phase_eigenvectors : numpy.ndarray
        Filtered phase-mode eigenvectors.
    glimmering_amplitude_eigenvalues : numpy.ndarray
        Weaker amplitude-mode eigenvalues from near-degenerate doublets.
    glimmering_amplitude_weights : numpy.ndarray
        Weaker amplitude-mode weights from near-degenerate doublets.
    glimmering_amplitude_eigenvectors : numpy.ndarray
        Weaker amplitude-mode eigenvectors from near-degenerate doublets.
    glimmering_phase_eigenvalues : numpy.ndarray
        Weaker phase-mode eigenvalues from near-degenerate doublets.
    glimmering_phase_weights : numpy.ndarray
        Weaker phase-mode weights from near-degenerate doublets.
    glimmering_phase_eigenvectors : numpy.ndarray
        Weaker phase-mode eigenvectors from near-degenerate doublets.

    Notes
    -----
    Eigenvectors are sign-normalized after filtering. Amplitude eigenvectors
    are flipped if the sum over their first ``N`` components is negative, while
    phase eigenvectors are flipped if their total sum is negative.
    """

    def __init__(self, system_data, epsilon_space):
        """Initialize the purger and process full-diagonalization data.

        Parameters
        ----------
        system_data : mapping
            Dictionary-like object containing raw amplitude and phase
            eigenvalues, weights, eigenvectors, and continuum boundaries.
        epsilon_space : array-like
            Grid used for plotting and for splitting amplitude eigenvectors
            into two components of length ``N``.
        """
        temp_amplitude_evs, temp_amplitude_weights, temp_amplitude_vectors = _fill_temporaries(
            system_data["amplitude.eigenvalues"], 
            system_data["amplitude.weights"][0], 
            system_data["amplitude.first_eigenvectors"],
            system_data["continuum_boundaries"][0]
        )
        temp_phase_evs, temp_phase_weights, temp_phase_vectors = _fill_temporaries(
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
        """Plot a normalized eigenvector line.

        For vectors of length ``2 * N``, the first and second halves are plotted
        on two separate axes. For vectors of any other length, the vector is
        plotted on a single axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or sequence of matplotlib.axes.Axes
            Axis or axes on which to plot. If ``len(y) == 2 * N``, ``ax`` is
            expected to contain two axes. Otherwise, ``ax`` is expected to be a
            single axis.
        y : array-like
            Eigenvector or vector-like data to plot.
        combined_norm : bool, default=True
            If ``True`` and ``len(y) == 2 * N``, both halves are normalized by
            the same maximum absolute value. If ``False``, both halves are
            normalized separately.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        float
            Normalization factor used for plotting. If ``combined_norm`` is
            ``False`` for a two-component vector, the returned value is
            ``norm1 / norm2``.

        Notes
        -----
        The plotted data are normalized by the inverse maximum absolute value.
        If the input vector contains only zeros, this may result in division by
        zero.
        """
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
        """Shared implementation for plotting selected eigenvectors.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or sequence of matplotlib.axes.Axes
            Axis or axes on which to plot.
        eigenvectors : array-like
            Collection of eigenvectors.
        eigenvalues : array-like
            Eigenvalues associated with ``eigenvectors``.
        which : int or {'all'}
            Which eigenvector to plot. If ``'all'``, all available eigenvectors
            are plotted.
        label_energy : bool
            If ``True``, include eigenvalue information in the plot label.
        combined_norm : bool, default=True
            Normalization mode forwarded to :meth:`plot_line`.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        None
            The method modifies the provided axes in place.
        """
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
        """Plot filtered phase-mode eigenvectors.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis on which to plot.
        which : int or {'all'}, default='all'
            Which phase eigenvector to plot. If ``'all'``, all filtered phase
            eigenvectors are plotted.
        label_energy : bool, default=False
            If ``True``, include eigenvalue information in the plot label.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        None
            The method modifies ``ax`` in place.
        """
        self.__plot_impl__(ax, self.phase_eigenvectors, self.phase_eigenvalues, which, label_energy, **plot_kwargs)
    
    ##############################################################
    def plot_glimmering_phase(self, ax, which='all', label_energy=False, **plot_kwargs):
        """Plot glimmering phase-mode eigenvectors.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis on which to plot.
        which : int or {'all'}, default='all'
            Which glimmering phase eigenvector to plot. If ``'all'``, all
            glimmering phase eigenvectors are plotted.
        label_energy : bool, default=False
            If ``True``, include eigenvalue information in the plot label.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        None
            The method modifies ``ax`` in place.
        """
        self.__plot_impl__(ax, self.glimmering_phase_eigenvectors, self.glimmering_phase_eigenvalues, which, label_energy, **plot_kwargs)
    
    ##############################################################
    def plot_amplitude(self, axes, which='all', label_energy=False, combined_norm=True, **plot_kwargs):
        """Plot filtered amplitude-mode eigenvectors.

        Amplitude eigenvectors are usually two-component vectors of length
        ``2 * N`` and are therefore plotted on two axes.

        Parameters
        ----------
        axes : sequence of matplotlib.axes.Axes
            Axes on which to plot the two amplitude components.
        which : int or {'all'}, default='all'
            Which amplitude eigenvector to plot. If ``'all'``, all filtered
            amplitude eigenvectors are plotted.
        label_energy : bool, default=False
            If ``True``, include eigenvalue information in the plot label.
        combined_norm : bool, default=True
            If ``True``, both amplitude components use one common
            normalization. If ``False``, each component is normalized
            separately.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        None
            The method modifies ``axes`` in place.
        """
        self.__plot_impl__(axes, self.amplitude_eigenvectors, self.amplitude_eigenvalues, which, 
                           label_energy, combined_norm, **plot_kwargs)
    
    ##############################################################
    def plot_glimmering_amplitude(self, axes, which='all', label_energy=False, combined_norm=True, **plot_kwargs):
        """Plot glimmering amplitude-mode eigenvectors.

        Parameters
        ----------
        axes : sequence of matplotlib.axes.Axes
            Axes on which to plot the two amplitude components.
        which : int or {'all'}, default='all'
            Which glimmering amplitude eigenvector to plot. If ``'all'``, all
            glimmering amplitude eigenvectors are plotted.
        label_energy : bool, default=False
            If ``True``, include eigenvalue information in the plot label.
        combined_norm : bool, default=True
            If ``True``, both amplitude components use one common
            normalization. If ``False``, each component is normalized
            separately.
        **plot_kwargs
            Additional keyword arguments forwarded to Matplotlib ``plot``.

        Returns
        -------
        None
            The method modifies ``axes`` in place.
        """
        self.__plot_impl__(axes, self.glimmering_amplitude_eigenvectors, self.glimmering_amplitude_eigenvalues, which, 
                           label_energy, combined_norm, **plot_kwargs)
        
    # returns a tuple (sum_k alpha_k^2, sum_k nu_k^2)
    def integral_amplitude(self, which):
        """Compute squared-component integrals of an amplitude eigenvector.

        The selected amplitude eigenvector is split into two components of
        length ``N``. The method returns the sum of squared entries of each
        component.

        Parameters
        ----------
        which : int
            Index of the amplitude mode.

        Returns
        -------
        tuple[float, float]
            Tuple ``(sum_k alpha_k^2, sum_k nu_k^2)`` containing the squared
            norms of the first and second halves of the selected amplitude
            eigenvector.

        Raises
        ------
        ValueError
            If ``which`` does not refer to an available amplitude eigenvector.
        """
        if which < len(self.amplitude_eigenvectors):
            pc = np.sum(self.amplitude_eigenvectors[which][:self.N]**2)
            n  = np.sum(self.amplitude_eigenvectors[which][self.N:]**2)
            #print(f"Pair creation integral: {pc}\nOccupation integral: {n}")
            #print(f"Combined: {pc + n}")
            return pc, n
        else:
            raise ValueError(f"Requested integral for amplitude mode {which}, but only {len(self.amplitude_eigenvectors)} available.")
        
    def to_dict(self, suffix):
        """Export selected filtered eigenvalues and eigenvectors as a dictionary.

        Parameters
        ----------
        suffix : str
            Suffix appended to each output key.

        Returns
        -------
        dict
            Dictionary containing filtered amplitude and phase eigenvalues and
            eigenvectors with keys of the form:

            - ``f"energies_higgs_{suffix}"``
            - ``f"energies_phase_{suffix}"``
            - ``f"amplitudes_higgs_{suffix}"``
            - ``f"amplitudes_phase_{suffix}"``
        """
        return { 
                f"energies_higgs_{suffix}" : self.amplitude_eigenvalues,
                f"energies_phase_{suffix}" : self.phase_eigenvalues,
                f"amplitudes_higgs_{suffix}" : self.amplitude_eigenvectors,
                f"amplitudes_phase_{suffix}" : self.phase_eigenvectors
            }