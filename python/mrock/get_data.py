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
    return (float(obj) if not isinstance(obj, str) else obj)

def _convert_to_numpy(obj):
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: _convert_to_numpy(v) for k, v in obj.items()}
    return obj

class DataLoader:
    def __init__(self, data_dir=None):
        if data_dir is None:
            self.data_dir = os.environ.get("MROCK_DATA_DIR", DEFAULT_DATA_DIR)
        else:    
            self.data_dir = Path(data_dir).expanduser().resolve()

    def _to_path(self, model, subdir, **kwargs):
        base = self.data_dir / model / Path(subdir)

        if kwargs:
            params = Path(*[f"{key}={value}" for key, value in kwargs.items()])
            return base / params

        return base

    def load_panda_file(self, file, numpy_conversion=True):
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
        path = self._to_path(model, subdir, **kwargs) / file
        data = self.load_panda_file(path, numpy_conversion=numpy_conversion)

        if print_date:
            print(f"Loaded data has been produced on {data['time']}")

        return data

    def get_all_files(self, folder, pattern):
        folder = Path(folder)
        search_obj = self.data_dir / folder / "**" / pattern
        return glob.glob(str(search_obj), recursive=True)

    def load_all(self, folder, pattern, condition=None):
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
        path = self._to_path(model, subdir) / file
        return pd.read_pickle(path)

def continuum_params(N_k, T, coulomb_scaling, screening, k_F, g, omega_D):
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
    return  {
                "T": _float_or_str(T), 
                "U": _float_or_str(U), 
                "V": _float_or_str(V)
            }

def hhg_params(T, E_F, v_F, band_width, field_amplitude, photon_energy, tau_diag, tau_offdiag, t0):
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
    return  {
                "N" : N,
                "g" : _float_or_str(g),
                "U" : _float_or_str(U),
                "E_F" : _float_or_str(E_F),
                "omega_D" : _float_or_str (omega_D)
            }

def dwave_params(N, g, V, E_F, omega_D):
    return  {
                "N" : N,
                "V" : _float_or_str(V),
                "E_F" : _float_or_str(E_F),
                "g" : _float_or_str(g),
                "omega_D" : _float_or_str (omega_D)
            }