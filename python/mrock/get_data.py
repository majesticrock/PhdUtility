import numpy as np
import gzip
import json
import pandas as pd
import glob
import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def __to_path(model, subfolder, **kwargs):
    if len(kwargs) != 0:
        params = os.path.join(*(f"{key}={value}" for key, value in kwargs.items()))
        return os.path.join(CURRENT_DIR, model, os.path.normpath(subfolder), params)
    return os.path.join(CURRENT_DIR, model, os.path.normpath(subfolder))

def convert_to_numpy(obj):
    if isinstance(obj, list):
        return np.array(obj)
    if isinstance(obj, dict):
        return {k: convert_to_numpy(v) for k, v in obj.items()}
    return obj

def load_panda_file(file, numpy_conversion=True):
    with gzip.open(file, 'rt') as f_open:
        if numpy_conversion:
            jData = json.load(f_open, object_hook=convert_to_numpy)
        else:
            jData = json.load(f_open)
    
    if 'data' in jData:
        main_data = {k: v for k, v in jData.items() if k != 'data'}
        main_df = pd.Series(main_data)
        main_df['data'] = pd.DataFrame(jData['data'])
    else:
        main_df = pd.json_normalize(jData, max_level=1).iloc[0]
    return main_df

def load_panda(model, subfolder, file, print_date=True, numpy_conversion=True, **kwargs):
    data = load_panda_file(os.path.join(__to_path(model, subfolder, **kwargs), file), numpy_conversion=numpy_conversion)
    if print_date:
        print(f"Loaded data has been produced on {data['time']}")
    return data

def __float_or_str__(obj):
    return (float(obj) if not isinstance(obj, str) else obj)

def continuum_params(N_k, T, coulomb_scaling, screening, k_F, g, omega_D):
    cp_screening = 0.0 if abs(coulomb_scaling) < 1e-12 else screening
    return  {
                "N_k"             : int(N_k) if not isinstance(N_k, str) else N_k, 
                "T"               : __float_or_str__(T),
                "coulomb_scaling" : __float_or_str__(coulomb_scaling),
                "screening"       : __float_or_str__(cp_screening),
                "k_F"             : __float_or_str__(k_F),
                "g"               : __float_or_str__(g), 
                "omega_D"         : __float_or_str__(omega_D)
            }

def hubbard_params(T, U, V):
    return  {
                "T": __float_or_str__(T), 
                "U": __float_or_str__(U), 
                "V": __float_or_str__(V)
            }

def hhg_params(T, E_F, v_F, band_width, field_amplitude, photon_energy, tau_diag, tau_offdiag, t0):
    return  {
                "T": __float_or_str__(T),
                "E_F": __float_or_str__(E_F),
                "v_F": __float_or_str__(v_F),
                "band_width": __float_or_str__(band_width),
                "field_amplitude": __float_or_str__(field_amplitude),
                "photon_energy": __float_or_str__(photon_energy),
                "tau_diag": __float_or_str__(tau_diag),
                "tau_offdiag" : __float_or_str__(tau_offdiag),
                "t0" : __float_or_str__(t0)
            }

def lattice_cut_params(N, g, U, E_F, omega_D):
    return  {
                "N" : N,
                "g" : __float_or_str__(g),
                "U" : __float_or_str__(U),
                "E_F" : __float_or_str__(E_F),
                "omega_D" : __float_or_str__ (omega_D)
            }

def dwave_params(N, g, V, E_F, omega_D):
    return  {
                "N" : N,
                "V" : __float_or_str__(V),
                "E_F" : __float_or_str__(E_F),
                "g" : __float_or_str__(g),
                "omega_D" : __float_or_str__ (omega_D)
            }

def get_all_files(folder, pattern):
    if isinstance(folder, str):
        folder = os.path.normpath(folder)
    search_obj = os.path.join(CURRENT_DIR, folder, os.path.normpath("**/"), pattern)
    return glob.glob(search_obj, recursive=True)

def load_all(folder, pattern, condition=None):
    all_files = get_all_files(folder, pattern)
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
    return pd.concat([load_panda_file(file) for file in all_files], ignore_index=True, axis=1).T.convert_dtypes()

def load_pickle(folder, file):
    if isinstance(folder, str):
        folder = os.path.normpath(folder)
    return pd.read_pickle(os.path.join(CURRENT_DIR, folder, file))

#pd_data = load_all("continuum/disc_2000", "gap.json.gz")
#print(pd_data["data"][0]["Delta_Coulomb"] + pd_data["data"][0]["Delta_Phonon"])