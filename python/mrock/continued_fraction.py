"""Continued-fraction utilities for resolvent and spectral-density calculations.

This module provides a high-level Python wrapper around the compiled
``cpp_continued_fraction`` backend. It loads continued-fraction coefficients
from a pandas-like data object, determines a suitable termination depth, and
evaluates quantities such as resolvents, spectral densities, denominators, and
bound-state classifications.

The main entry point is :class:`ContinuedFraction`.
"""

import numpy as np
import matplotlib.pyplot as plt

from . import cpp_continued_fraction as ccf

NORM_FACTOR = -(1. / np.pi)


class ContinuedFraction:
    """Evaluate continued fractions from stored resolvent data.

    This class wraps continued-fraction coefficient data contained in a
    pandas-like data structure and delegates numerical evaluations to the
    compiled ``cpp_continued_fraction`` backend.

    Parameters
    ----------
    data_frame : pandas.Series or pandas.DataFrame-like
        Data object containing continued-fraction data. It is expected to
        contain a ``"continuum_boundaries"`` entry with exactly two values and
        entries of the form ``"resolvents.<resolvent_name>"``.
        
        If you generated your data with the iEoM library and used the get_data script
        with DataLoader.load_panda(...) you will get an object that you can
        pass straight to this class.
    messages : bool, default=True
        If ``True``, print diagnostic messages such as the chosen termination
        depth.
        
    ignore_first : int, default=5
        Number of initial coefficients to ignore when searching for the best
        termination depth.
    ignore_last : int, default=80
        Maximum number of coefficients from the beginning of the coefficient
        sequence to inspect when estimating the termination depth.
        
        There is no general way to choose ignore_first and ignore_last.
        If the continued fraction terminates too early, you will miss out on
        physical features such as bound states and collective modes.
        However, due to numerical constraints, the coefficients accumulate
        errors the deeper you go, and you might see features akin to overfitting.
        Generally, it seems that the larger the system (or discretization) is the
        depper you can and should go with the continued fraction.
        Nonetheless, the proper depth will depend on the studied system.

    Attributes
    ----------
    messages : bool
        Whether diagnostic messages are printed.
    ignore_first : int
        Number of initial coefficients ignored during termination-depth search.
    ignore_last : int
        Maximum number of coefficients inspected during termination-depth
        search.
    continuum_boundaries_squared : numpy.ndarray
        Squared lower and upper continuum boundaries.
    a_infinity : float
        Asymptotic value of the continued-fraction ``a`` coefficients.
    b_infinity : float
        Asymptotic value derived from the continuum boundaries.
    data : pandas.Series or pandas.DataFrame-like
        Input data object containing continued-fraction coefficients.
    termination_resolvent_index : dict
        Cache mapping resolvent names and indices to selected termination
        depths.

    Raises
    ------
    ValueError
        If ``continuum_boundaries`` does not contain exactly two values or if
        the lower boundary is not smaller than the upper boundary.
    """

    def __init__(self, data_frame, messages=True, ignore_first=5, ignore_last=80):
        """Initialize a continued-fraction evaluator.
        See the class description for more details.
        
        Parameters
        ----------
        data_frame : pandas.Series or pandas.DataFrame-like
            Data object containing the continuum boundaries and resolvent
            coefficient data.
        messages : bool, default=True
            If ``True``, print diagnostic messages.
        ignore_first : int, default=5
            Number of initial coefficients ignored when selecting the
            termination depth.
        ignore_last : int, default=80
            Maximum number of coefficients inspected when selecting the
            termination depth.
        """
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

    def terminator(self, omega, resolvent_name, resolvent_index=0, with_terminator=True):
        """Evaluate the terminator contribution for a resolvent.

        Parameters
        ----------
        omega : scalar or array-like
            Frequency or complex frequency values at which to evaluate the
            terminator.
        resolvent_name : str
            Name of the resolvent inside the data object.
        resolvent_index : int, default=0
            Index of the resolvent entry to use.
        with_terminator : bool, default=True
            Whether the returned continued-fraction data should include the
            terminator.

        Returns
        -------
        numpy.ndarray
            Terminator values evaluated at ``omega``.
        """
        data = self._get_cf_data(resolvent_name, resolvent_index, with_terminator)
        omega = np.asarray(omega, dtype=np.float64)
        return ccf.terminator(omega, data)
    
    def _coeffs_A(self, resolvent_name, resolvent_index):
        """Return the ``a_i`` continued-fraction coefficients.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.

        Returns
        -------
        numpy.ndarray
            Array containing the ``a_i`` coefficients.
        """
        return np.asarray(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["a_i"])
    
    def _coeffs_B(self, resolvent_name, resolvent_index):
        """Return the ``b_i`` continued-fraction coefficients.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.

        Returns
        -------
        numpy.ndarray
            Array containing the ``b_i`` coefficients.
        """
        return np.asarray(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["b_i"])
    
    def _len_A(self, resolvent_name, resolvent_index):
        """Return the number of ``a_i`` coefficients.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.

        Returns
        -------
        int
            Number of available ``a_i`` coefficients.
        """
        return len(self.data[f"resolvents.{resolvent_name}"][resolvent_index]["a_i"])
    
    def find_termination_depth(self, resolvent_name, resolvent_index):
        """Determine and cache a termination depth for a resolvent.

        The termination depth is chosen by comparing the continued-fraction
        coefficients against their asymptotic values and selecting the index
        with the smallest relative deviation after skipping the first
        ``ignore_first`` entries.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.

        Raises
        ------
        ValueError
            If ``ignore_first`` is greater than or equal to the effective
            ``ignore_last`` value.
        """
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
        """Ensure that a termination depth exists for a resolvent.

        If no cached termination-depth array exists for ``resolvent_name``, it
        is initialized. If the selected ``resolvent_index`` has no cached
        termination depth yet, :meth:`find_termination_depth` is called.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.
        """
        if not resolvent_name in self.termination_resolvent_index:
            self.termination_resolvent_index[resolvent_name] = np.zeros(len(self.data[f"resolvents.{resolvent_name}"]), dtype=int)
        if self.termination_resolvent_index[resolvent_name][resolvent_index] == 0:
            self.find_termination_depth(resolvent_name, resolvent_index)
    
    def _get_cf_data(self, resolvent_name, resolvent_index, with_terminator):
        """Construct backend continued-fraction data.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int
            Index of the resolvent entry.
        with_terminator : bool
            Whether to include the terminator in the backend data object.

        Returns
        -------
        cpp_continued_fraction.ContinuedFractionData
            Data object passed to the compiled continued-fraction backend.
        """
        self._termination_depth_if_required(resolvent_name, resolvent_index)
        return ccf.ContinuedFractionData(self.a_infinity, self.b_infinity**2, np.asarray(self.continuum_boundaries_squared), 
                                         np.asarray(self._coeffs_A(resolvent_name, resolvent_index)), np.asarray(self._coeffs_B(resolvent_name, resolvent_index)), 
                                         self.termination_resolvent_index[resolvent_name][resolvent_index], with_terminator)
    
    def continued_fraction(self, omega, resolvent_name, resolvent_index=0, with_terminator=True):
        """Evaluate a continued fraction based on the coefficients provided when constructing 
        an instance of this class.

        Parameters
        ----------
        omega : scalar or array-like
            Complex frequency value or values at which to evaluate the
            continued fraction.
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        complex or numpy.ndarray
            Continued-fraction value. A scalar is returned if ``omega`` is a
            scalar; otherwise an array is returned.
        """
        scalar_input = np.isscalar(omega)
        omega = np.atleast_1d(np.asarray(omega, dtype=np.complex128))
        result = ccf.continued_fraction(omega, self._get_cf_data(resolvent_name, resolvent_index, with_terminator))
        return result[0] if scalar_input else result
    
    def continued_fraction_varied_depth(self, omega, resolvent_name, shift_range, resolvent_index=0, with_terminator=True):
        """Evaluate the many continued fraction with varied termination depths.

        Parameters
        ----------
        omega : scalar or array-like
            Complex frequency value or values at which to evaluate the
            continued fraction.
        resolvent_name : str
            Name of the resolvent.
        shift_range : array-like
            Range of shifts applied to the selected termination depth.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        complex or numpy.ndarray
            Continued-fraction values computed with varied termination depths.
            A scalar-like first entry is returned if ``omega`` is a scalar;
            otherwise an array is returned.
        """
        scalar_input = np.isscalar(omega)
        omega = np.atleast_1d(np.asarray(omega, dtype=np.complex128))
        result = ccf.continued_fraction_varied_depth(omega, self._get_cf_data(resolvent_name, resolvent_index, with_terminator), shift_range)
        return result[0] if scalar_input else result

    def spectral_density(self, omega, resolvent_name, resolvent_index=0, with_terminator=True):
        """Evaluate the spectral density.

        The spectral density is computed from the imaginary part of the
        continued fraction using the normalization factor $1/\pi$.

        Parameters
        ----------
        omega : scalar or array-like
            Frequency or complex frequency value or values.
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        float or numpy.ndarray
            Spectral-density value or values.
        """
        return NORM_FACTOR * self.continued_fraction(omega, resolvent_name, resolvent_index, with_terminator).imag
    
    def spectral_density_varied_depth(self, omega, resolvent_name, shift_range, resolvent_index=0, with_terminator=True):
        """Evaluate many spectral densities with varied termination depths.

        Parameters
        ----------
        omega : scalar or array-like
            Frequency or complex frequency value or values.
        resolvent_name : str
            Name of the resolvent.
        shift_range : array-like
            Range of shifts applied to the selected termination depth.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        float or numpy.ndarray
            Spectral-density values computed with varied termination depths.
        """
        return NORM_FACTOR * self.continued_fraction_varied_depth(omega, resolvent_name, shift_range, resolvent_index, with_terminator).imag
    
    def real_part(self, omega, resolvent_name, resolvent_index=0, with_terminator=True):
        """Evaluate the real part of the continued fraction.

        Parameters
        ----------
        omega : scalar or array-like
            Complex frequency value or values.
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        float or numpy.ndarray
            Real part of the continued-fraction value.
        """
        return self.continued_fraction(omega, resolvent_name, resolvent_index, with_terminator).real
    
    def denominator(self, omega, resolvent_name, resolvent_index=0, with_terminator=True):
        """Evaluate the denominator of the continued fraction.

        Parameters
        ----------
        omega : scalar or array-like
            Complex frequency value or values.
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        with_terminator : bool, default=True
            Whether to use the asymptotic terminator.

        Returns
        -------
        complex or numpy.ndarray
            Denominator value. A scalar is returned if ``omega`` is a scalar;
            otherwise an array is returned.
        """
        scalar_input = np.isscalar(omega)
        omega = np.atleast_1d(np.asarray(omega, dtype=np.complex128))
        result = ccf.denominator(omega, self._get_cf_data(resolvent_name, resolvent_index, with_terminator))
        return result[0] if scalar_input else result
    
    # Does not really work for Goldstone peaks
    def classify_bound_states(self, resolvent_name, resolvent_index=0, n_scan=1000, weight_domega=1e-8, tolerance_bits=48, max_iter=200, is_phase_peak=None):
        """Classify bound states associated with a resolvent.

        This method delegates the bound-state search and classification to the
        compiled backend. Optionally, detected phase peaks can be filtered out.
        
        Atleast in the cases I considered, this method has trouble with Goldstone modes.
        These appear as delta'-peaks in the spectral function, i.e., as double-poles
        of the nature 1/omega^2 in the complex Green's function.
        In contrast, finite-energy subgap peaks are simple poles, for which the algorithm
        simply computes the residue. This works splendidly.

        Parameters
        ----------
        resolvent_name : str
            Name of the resolvent.
        resolvent_index : int, default=0
            Index of the resolvent entry.
        n_scan : int, default=1000
            Number of scan points used by the backend search.
        weight_domega : float, default=1e-8
            Small frequency offset used for estimating bound-state weights.
        tolerance_bits : int, default=48
            Precision target used by the backend root-finding routine.
        max_iter : int, default=200
            Maximum number of iterations used by the backend.
        is_phase_peak : callable, optional
            Predicate used to identify phase peaks. If provided, items for
            which ``is_phase_peak(item[0])`` returns ``True`` are removed from
            the result.

        Returns
        -------
        list
            Bound-state classification results returned by the backend,
            optionally filtered.

        Notes
        -----
        The existing implementation comment indicates that this method does not
        work reliably for Goldstone peaks.
        """
        result = ccf.classify_bound_states(self._get_cf_data(resolvent_name, resolvent_index, with_terminator=True), n_scan, weight_domega, tolerance_bits, max_iter)
        if is_phase_peak is not None:
            result = [ item for item in result if not is_phase_peak(item[0]) ]
        return result
                
    
    def mark_continuum(self, ax=None, scale_factor=1., label="Continuum"):
        """Mark the continuum interval in a Matplotlib plot.

        The continuum region is drawn as a shaded vertical span between the
        lower and upper continuum edges.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object on which to draw the continuum region. If ``None``, the
            current Matplotlib axes are used via :func:`matplotlib.pyplot.axvspan`.
        scale_factor : float, default=1.0
            Factor by which the continuum edges are multiplied before plotting.
        label : str or None, default="Continuum"
            Label assigned to the shaded region. If ``None``, no label is
            passed to Matplotlib.

        Returns
        -------
        None
            This method modifies the plot in place.
        """
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
        """Return the lower and upper continuum edges.

        Returns
        -------
        numpy.ndarray
            Array containing the lower and upper continuum boundaries. These
            are obtained by taking the square root of
            ``continuum_boundaries_squared``.
        """
        return np.sqrt(self.continuum_boundaries_squared)