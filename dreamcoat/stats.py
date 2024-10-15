import warnings
import numpy as np
from scipy.special import factorial


def std_bias_correction(sample_size):
    """Find bias correction for standard deviations computed from small numbers of
    normally distributed measurements.

    Parameters
    ----------
    sample_size : int
        Number of measurements from which the standard deviation was computed.

    Returns
    -------
    c4 : float
        The correction factor c4, based on the equations from
        https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
        The standard deviation obtained from np.std(a, ddof=1) should be divided by
        this c4 value to find the unbiased estimate.
    """
    k = np.floor(sample_size / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        c4 = np.where(
            sample_size % 2 == 0,
            # For even sample_size:
            np.sqrt(2 / (np.pi * (2 * k - 1)))
            * (2 ** (2 * k - 2) * factorial(k - 1) ** 2)
            / factorial(2 * k - 2),
            # For odd sample_size:
            np.sqrt(np.pi / k)
            * factorial(2 * k - 1)
            / (2 ** (2 * k - 1) * factorial(k - 1) ** 2),
        )
    return c4


def std_unbiased(a, axis=None, **kwargs):
    """Compute standard deviation unbiased for the number of samples, assuming a normal
    distribution.  The correction factor c4 is used, based on the equations from
    https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation

    Parameters
    ----------
    a : array-like
        The values to find the standard deviation of.
    axis : int, optional
        Which axis to find the standard deviation along, by default None.  For example,
        with an array of shape (1000, 4) and axis=1, this represents 1000 different
        groups of 4 measurements.

    Returns
    -------
    float
        The unbiased standard deviation, corrected for the number of samples.
    """
    if axis is None:
        size_a = np.size(a)
    else:
        size_a = np.shape(a)[axis]
    c4 = std_bias_correction(size_a)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice.")
        std_biased = np.nanstd(a, axis=axis, ddof=1, **kwargs)
        # ^ ddof must be 1 in order for the c4 correction to be valid
    return std_biased / c4
