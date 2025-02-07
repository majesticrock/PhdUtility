import numpy as np
from uncertainties import ufloat

def fit_result_to_ufloat(popt, pcov, i=None):
    if i is None:
        return [ufloat(popt[j], np.sqrt(pcov[j][j])) for j in range(len(popt))]
    return ufloat(popt[i], np.sqrt(pcov[i][i]))

def fit_result_to_string(popt, pcov, i):
    __u = fit_result_to_ufloat(popt, pcov, i)
    return '{:L}'.format(__u)