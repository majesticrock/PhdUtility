"""Utilities for loading compressed JSON and pickle data files.

This module provides a :class:`DataLoader` for locating and loading data files
from a configurable data directory. It also contains helper functions for
constructing parameter dictionaries used to identify or organize simulation
outputs.

The default data directory is resolved relative to this file, but it can be
overridden with the ``MROCK_DATA_DIR`` environment variable or by passing an
explicit path to :class:`DataLoader`.
"""

import numpy as np
import gzip
import json
import pandas as pd
import glob
import os
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = CURRENT_DIR.parents[2] / "data"


def _float_or_str(obj):
    """Convert an object to ``float`` unless it is already a string.

    Parameters
    ----------
    obj : Any
        Object to convert.

    Returns
    -------
    float or str
        ``float(obj)`` if ``obj`` is not a string, otherwise ``obj`` unchanged.

    Notes
    -----
    This helper is useful for building parameter dictionaries where numeric
    values should be normalized to floats, while special string values should
    be preserved.
    """
    return (float(obj) if not isinstance(obj, str) else obj)


def _convert_to_numpy(obj):
    """Recursively convert lists inside an object to NumPy arrays.

    Parameters
    ----------
    obj : Any
        Object returned while decoding JSON data. Lists are converted to
        :class:`numpy.ndarray` objects, dictionaries are processed recursively,
        and all other objects are returned unchanged.

    Returns
    -------
    Any
        The converted object.

    Examples
    --------
    ``[1, 2, 3]`` becomes ``np.array([1, 2, 3])``.

    ``{"x": [1, 2]}`` becomes ``{"x": np.array([1, 2])}``.
    """
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: _convert_to_numpy(v) for k, v in obj.items()}
    return obj


class DataLoader:
    """Loader for simulation data stored below a configurable data directory.

    The loader supports compressed JSON files, which are converted into pandas
    objects, and pickle files loaded through :func:`pandas.read_pickle`.

    Parameters
    ----------
    data_dir : str or pathlib.Path, optional
        Base directory containing the data. If omitted, the environment variable
        ``MROCK_DATA_DIR`` is used when available. Otherwise,
        :data:`DEFAULT_DATA_DIR` is used.

    Attributes
    ----------
    data_dir : pathlib.Path or str
        Base directory used to resolve data file paths.
    """

    def __init__(self, data_dir=None):
        """Initialize a :class:`DataLoader`.

        Parameters
        ----------
        data_dir : str or pathlib.Path, optional
            Base data directory. If ``None``, the value of the
            ``MROCK_DATA_DIR`` environment variable is used if present;
            otherwise :data:`DEFAULT_DATA_DIR` is used.
        """
        if data_dir is None:
            self.data_dir = os.environ.get("MROCK_DATA_DIR", DEFAULT_DATA_DIR)
        else:
            self.data_dir = Path(data_dir).expanduser().resolve()

    def _to_path(self, model, subdir, **kwargs):
        """Construct a path below the data directory.

        Parameters
        ----------
        model : str or pathlib.Path
            Name of the model directory below ``data_dir``.
        subdir : str or pathlib.Path
            Subdirectory below the model directory.
        **kwargs
            Optional key-value pairs appended as nested path components in the
            form ``key=value``. The order follows the insertion order of the
            provided keyword arguments.

        Returns
        -------
        pathlib.Path
            Constructed path.

        Examples
        --------
        ``_to_path("hubbard", "runs", T=0.1, U=2.0)`` returns a path like::

            data_dir / "hubbard" / "runs" / "T=0.1" / "U=2.0"
        """
        base = self.data_dir / model / Path(subdir)

        if kwargs:
            params = Path(*[f"{key}={value}" for key, value in kwargs.items()])
            return base / params

        return base

    def load_panda_file(self, file, numpy_conversion=True):
        """Load a compressed JSON file into a pandas object.

        Parameters
        ----------
        file : str or pathlib.Path
            Path to a gzip-compressed JSON file.
        numpy_conversion : bool, default=True
            If ``True``, lists in the decoded JSON data are recursively
            converted to NumPy arrays using :func:`_convert_to_numpy`.

        Returns
        -------
        pandas.Series
            A pandas Series containing the loaded metadata and data. If the JSON
            object contains a ``"data"`` key, all other fields are stored as
            scalar entries and ``"data"`` is converted to a
            :class:`pandas.DataFrame`. Otherwise, the JSON object is normalized
            using :func:`pandas.json_normalize`.

        Notes
        -----
        The input file is expected to be gzip-compressed and readable in text
        mode.
        """
        with gzip.open(file, "rt") as f_open:
            if numpy_conversion:
                jData = json.load(f_open, object_hook=_convert_to_numpy)
            else:
                jData = json.load(f_open)

        if "data" in jData:
            main_data = {k: v for k, v in jData.items() if k != "data"}
            main_df = pd.Series(main_data)
            main_df["data"] = pd.DataFrame(jData["data"])
        else:
            main_df = pd.json_normalize(jData, max_level=1).iloc[0]

        return main_df

    def load_panda(self, model, subdir, file, print_date=True, numpy_conversion=True, **kwargs):
        """Load a compressed JSON file from a model-specific data directory.

        Parameters
        ----------
        model : str or pathlib.Path
            Name of the model directory below ``data_dir``.
        subdir : str or pathlib.Path
            Subdirectory below the model directory.
        file : str or pathlib.Path
            File name to load.
        print_date : bool, default=True
            If ``True``, print the value stored under the ``"time"`` field of
            the loaded data.
        numpy_conversion : bool, default=True
            If ``True``, recursively convert JSON lists to NumPy arrays.
        **kwargs
            Optional path parameters appended as ``key=value`` subdirectories.

        Returns
        -------
        pandas.Series
            Loaded data object returned by :meth:`load_panda_file`.

        Raises
        ------
        KeyError
            If ``print_date`` is ``True`` and the loaded data does not contain
            a ``"time"`` field.
        """
        path = self._to_path(model, subdir, **kwargs) / file
        data = self.load_panda_file(path, numpy_conversion=numpy_conversion)

        if print_date:
            print(f"Loaded data has been produced on {data['time']}")

        return data

    def get_all_files(self, folder, pattern):
        """Find all files matching a pattern below a folder.

        Parameters
        ----------
        folder : str or pathlib.Path
            Folder below ``data_dir`` in which the recursive search starts.
        pattern : str
            Glob pattern used to match file names, for example ``"*.json.gz"``.

        Returns
        -------
        list[str]
            List of matching file paths as strings.

        Notes
        -----
        The search is recursive and uses the pattern::

            data_dir / folder / "**" / pattern
        """
        folder = Path(folder)
        search_obj = self.data_dir / folder / "**" / pattern
        return glob.glob(str(search_obj), recursive=True)

    def load_all(self, folder, pattern, condition=None):
        """Load and concatenate all matching compressed JSON files.

        Parameters
        ----------
        folder : str or pathlib.Path
            Folder below ``data_dir`` in which to search recursively.
        pattern : str
            Glob pattern used to match file names.
        condition : str or iterable of str, optional
            Optional directory-name filter. If provided, only files whose path
            contains the condition as a full path component are loaded. If an
            iterable is provided, all conditions are applied successively.

        Returns
        -------
        pandas.DataFrame
            DataFrame created by loading all matching files with
            :meth:`load_panda_file`, concatenating them column-wise, transposing
            the result, and applying :meth:`pandas.DataFrame.convert_dtypes`.

        Notes
        -----
        Path separators are handled separately for Windows and POSIX systems.
        """
        all_files = self.get_all_files(folder, pattern)
        if condition is not None:
            if not isinstance(condition, str):
                for cond in condition:
                    if os.name == 'nt':
                        cond = f"\\{cond}\\"
                    else:
                        cond = f"/{cond}/"
                    all_files = [file for file in all_files if cond in file]
            else:
                if os.name == 'nt':
                    cond = f"\\{condition}\\"
                else:
                    cond = f"/{condition}/"
                all_files = [file for file in all_files if cond in file]
        return pd.concat([self.load_panda_file(file) for file in all_files], ignore_index=True, axis=1).T.convert_dtypes()

    def load_pickle(self, model, subdir, file):
        """Load a pickle file from a model-specific data directory.

        Parameters
        ----------
        model : str or pathlib.Path
            Name of the model directory below ``data_dir``.
        subdir : str or pathlib.Path
            Subdirectory below the model directory.
        file : str or pathlib.Path
            Pickle file name to load.

        Returns
        -------
        Any
            Object loaded by :func:`pandas.read_pickle`.
        """
        path = self._to_path(model, subdir) / file
        return pd.read_pickle(path)


def continuum_params(N_k, T, coulomb_scaling, screening, k_F, g, omega_D):
    """Create a parameter dictionary for a continuum model calculation.

    Parameters
    ----------
    N_k : int or str
        Number of momentum points. Converted to ``int`` unless it is a string.
    T : float-like or str
        Temperature.
    coulomb_scaling : float-like or str
        Coulomb interaction scaling factor.
    screening : float-like or str
        Screening parameter. Set to ``0.0`` if ``coulomb_scaling`` is
        numerically close to zero.
    k_F : float-like or str
        Fermi momentum.
    g : float-like or str
        Phonon-mediated coupling constant.
    omega_D : float-like or str
        Debye frequency.

    Returns
    -------
    dict
        Dictionary containing normalized continuum-model parameters.

    Notes
    -----
    Numeric values are converted to ``float`` using :func:`_float_or_str`,
    except for ``N_k``, which is converted to ``int`` unless it is a string.
    """
    cp_screening = 0.0 if abs(coulomb_scaling) < 1e-12 else screening
    return  {
                "N_k"             : int(N_k) if not isinstance(N_k, str) else N_k, 
                "T"               : _float_or_str(T),
                "coulomb_scaling" : _float_or_str(coulomb_scaling),
                "screening"       : _float_or_str(cp_screening),
                "k_F"             : _float_or_str(k_F),
                "g"               : _float_or_str(g), 
                "omega_D"         : _float_or_str(omega_D)
            }


def hubbard_params(T, U, V):
    """Create a parameter dictionary for a Hubbard model calculation.

    Parameters
    ----------
    T : float-like or str
        Temperature.
    U : float-like or str
        On-site interaction strength.
    V : float-like or str
        Non-local interaction strength.

    Returns
    -------
    dict
        Dictionary containing normalized Hubbard-model parameters.
    """
    return  {
                "T": _float_or_str(T), 
                "U": _float_or_str(U), 
                "V": _float_or_str(V)
            }


def hhg_params(T, E_F, v_F, band_width, field_amplitude, photon_energy, tau_diag, tau_offdiag, t0):
    """Create a parameter dictionary for a high-harmonic-generation calculation.

    Parameters
    ----------
    T : float-like or str
        Temperature.
    E_F : float-like or str
        Fermi energy.
    v_F : float-like or str
        Fermi velocity.
    band_width : float-like or str
        Band width.
    field_amplitude : float-like or str
        Driving-field amplitude.
    photon_energy : float-like or str
        Photon energy.
    tau_diag : float-like or str
        Relaxation time for diagonal density-matrix components.
    tau_offdiag : float-like or str
        Relaxation time for off-diagonal density-matrix components.
    t0 : float-like or str
        Time delay between the two laser pulses.

    Returns
    -------
    dict
        Dictionary containing normalized high-harmonic-generation parameters.
    """
    return  {
                "T": _float_or_str(T),
                "E_F": _float_or_str(E_F),
                "v_F": _float_or_str(v_F),
                "band_width": _float_or_str(band_width),
                "field_amplitude": _float_or_str(field_amplitude),
                "photon_energy": _float_or_str(photon_energy),
                "tau_diag": _float_or_str(tau_diag),
                "tau_offdiag" : _float_or_str(tau_offdiag),
                "t0" : _float_or_str(t0)
            }


def lattice_cut_params(N, g, U, E_F, omega_D):
    """Create a parameter dictionary for a lattice-cut calculation.

    Parameters
    ----------
    N : int or str
        System size.
    g : float-like or str
        Phonon mediated coupling constant.
    U : float-like or str
        Hubbard-like interaction strength.
    E_F : float-like or str
        Fermi energy.
    omega_D : float-like or str
        Debye frequency.

    Returns
    -------
    dict
        Dictionary containing normalized lattice-cut parameters.
    """
    return  {
                "N" : N,
                "g" : _float_or_str(g),
                "U" : _float_or_str(U),
                "E_F" : _float_or_str(E_F),
                "omega_D" : _float_or_str (omega_D)
            }


def dwave_params(N, g, V, E_F, omega_D):
    """Create a parameter dictionary for a d-wave calculation.

    Parameters
    ----------
    N : int or str
        System size.
    g : float-like or str
        s-symmetry coupling constant.
    V : float-like or str
        d-symmetry coupling constant.
    E_F : float-like or str
        Fermi energy.
    omega_D : float-like or str
        Debye frequency.

    Returns
    -------
    dict
        Dictionary containing normalized d-wave-model parameters.
    """
    return  {
                "N" : N,
                "V" : _float_or_str(V),
                "E_F" : _float_or_str(E_F),
                "g" : _float_or_str(g),
                "omega_D" : _float_or_str (omega_D)
            }