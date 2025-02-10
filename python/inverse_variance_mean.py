import numpy as np
# Computes the mean of quantities with a given standard deviation. All standard deviations must be finite.
def inverse_variance_mean(nominal_values, standard_deviations):
    inv_variances = 1. / standard_deviations**2
    inv_sigma_sum = np.sum(inv_variances)
    weighted_sum = np.sum(nominal_values * inv_variances)
    return weighted_sum / inv_sigma_sum, 1. / np.sqrt(inv_sigma_sum)
    
    