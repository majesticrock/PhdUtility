import math
from uncertainties import ufloat

def fit_result_to_ufloat(popt, pcov, i=None):
    if i is None:
        return [ufloat(popt[j], math.sqrt(pcov[j][j])) for j in range(len(popt))]
    return ufloat(popt[i], math.sqrt(pcov[i][i]))

def format_with_uncertainty(x, DeltaX):
    DeltaX = abs(DeltaX)
    if DeltaX == 0:
        return f"{x:.10g} ± {DeltaX:.10g}"
    
    error_order = -int(math.floor(math.log10(DeltaX)))
    rounded_x = round(x, error_order)
    rounded_DeltaX = round(DeltaX, error_order)

    return f"{rounded_x:.{error_order}f} \\pm {rounded_DeltaX:.{error_order}f}"

def fit_result_to_string(popt, pcov, i):
    return format_with_uncertainty(popt[i], math.sqrt(pcov[i][i]))

if __name__ == '__main__':
    x = 123.45678
    DeltaX = 0.01234

    print(format_with_uncertainty(x, DeltaX))