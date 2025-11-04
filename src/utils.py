import numpy as np
from scipy.ndimage import gaussian_filter1d


def get_best_t(point, theory_xmins, theory_ymins, T_values, v0_values, x_weight=1e-4):
    """
    Find the best-fitting loop period T and effective fragment length v0
    for a given experimental minimum (x, y) in the (s_min, y_min) plane.

    Parameters
    ----------
    point : tuple or list
        Experimental minimum as (x, y), where:
        - x is s_min (in the same units as theory_xmins, typically kb)
        - y is y_min (log-derivative value at the minimum)
    theory_xmins : dict[int, np.ndarray]
        Dictionary mapping T -> array of s_min(T, v0) values
        evaluated on a grid of v0 values.
        Example: theory_xmins[150] has shape (len(v0_values),)
    theory_ymins : dict[int, np.ndarray]
        Dictionary mapping T -> array of y_min(T, v0) values
        evaluated on the same v0 grid as theory_xmins.
    T_values : array-like
        List or array of candidate loop periods T (in kb),
        e.g. np.arange(100, 301, 10).
    v0_values : array-like
        Grid of effective fragment lengths v0 (in kb) used to compute
        theory_xmins[T] and theory_ymins[T], e.g.
        [0.0001, 0.2, 0.4, ..., 14, 15].
    x_weight : float, optional
        Relative weight of deviations along the x-axis (s_min) in the
        error function. The error is defined as:
            err = x_weight * (x_th - x_exp)^2 + (y_th - y_exp)^2
        Default: 1e-4.

    Returns
    -------
    best_T : float
        Best-fitting loop period (in kb).
    best_v0 : float
        Best-fitting effective fragment length (in kb).
    best_err : float
        Minimum value of the weighted squared error.
    """
    x_exp, y_exp = float(point[0]), float(point[1])

    best_err = np.inf
    best_T = None
    best_v0 = None

    v0_values = np.asarray(v0_values)

    for T in T_values:
        x_theory = np.asarray(theory_xmins[T])
        y_theory = np.asarray(theory_ymins[T])

        # Weighted squared distance between (x_exp, y_exp)
        # and the theoretical curve at this T
        err_array = x_weight * (x_theory - x_exp) ** 2 + (y_theory - y_exp) ** 2
        idx = np.argmin(err_array)
        err_T = err_array[idx]

        if err_T < best_err:
            best_err = err_T
            best_T = T
            best_v0 = v0_values[idx]

    return best_T, best_v0, best_err


def get_theory_curves(
    T_values=None,
    v0_values=None,
    s_min=1.0,
    s_max=1e4,
    s_step=0.01,
    smooth_sigma=0.0,
):
    """
    Compute theoretical positions and depths of log-derivative minima
    for the short-range Hi-C signature, across a grid of loop periods T
    and effective fragment lengths v0.

    The theory is based on Eq. (3) in the main text:
        y(s) = d log P(s) / d log s
             = - 3γ / [2(γ + 1)] + (s / T) * 3γ(γ + 2) / [2(γ + 1)^2],
    where γ(s) = 2s / (3 v0).

    For each T and v0, we evaluate y(s) over s ∈ [s_min, s_max],
    optionally smooth y(s), and record:
        - s_min(T, v0)  = argmin_s y(s)
        - y_min(T, v0)  = min_s y(s)

    Parameters
    ----------
    T_values : array-like, optional
        Loop periods T (in kb) to evaluate. If None, defaults to
        np.arange(100, 301, 10) → 100, 110, ..., 300 kb.
    v0_values : array-like, optional
        Effective fragment lengths v0 (in kb) to evaluate. If None,
        uses the grid from the paper:
        [0.0001, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 2.5, 3, 3.5,
         4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.6,
         10, 10.5, 11, 11.5, 12, 12.5, 13, 14, 15]
    s_min : float, optional
        Minimum genomic separation (in kb) to evaluate y(s).
        Default: 1.0
    s_max : float, optional
        Maximum genomic separation (in kb) to evaluate y(s).
        Default: 1e4
    s_step : float, optional
        Step size for the s grid. Default: 0.01
    smooth_sigma : float, optional
        Standard deviation for Gaussian smoothing of y(s)
        (in units of index along s). If 0.0, no smoothing is applied.
        Default: 0.0

    Returns
    -------
    theory_xmins : dict[int, np.ndarray]
        Dictionary mapping T -> array of s_min(T, v0) for each v0 in v0_values.
        Shape of each value: (len(v0_values),)
    theory_ymins : dict[int, np.ndarray]
        Dictionary mapping T -> array of y_min(T, v0) for each v0 in v0_values.
        Shape of each value: (len(v0_values),)
    s_grid : np.ndarray
        The s grid (in kb) used for the calculations.
    v0_values : np.ndarray
        The v0 grid (in kb) actually used (as a NumPy array).
    """

    # Default grids if not provided
    if T_values is None:
        T_values = np.arange(100, 301, 10)  # 100, 110, ..., 300 kb
    if v0_values is None:
        v0_values = [
            0.0001, 0.2, 0.4, 0.6, 0.8, 1,
            1.5, 2, 2.5, 3, 3.5,
            4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5,
            8, 8.5, 9, 9.6,
            10, 10.5, 11, 11.5, 12, 12.5,
            13, 14, 15,
        ]

    T_values = np.asarray(T_values, dtype=float)
    v0_values = np.asarray(v0_values, dtype=float)

    # Genomic separation grid (in kb)
    s_grid = np.arange(s_min, s_max, s_step, dtype=float)

    # Containers for results
    theory_xmins = {}
    theory_ymins = {}

    # Loop over T values
    for T in T_values:
        s_min_list = []
        y_min_list = []

        # Loop over v0 values
        for v0 in v0_values:
            # γ(s) = 2s / (3 v0)
            gamma = 2.0 * s_grid / (3.0 * v0)

            # y(s) = -3γ / [2(γ+1)] + (s/T) * 3γ(γ+2) / [2(γ+1)^2]
            y = (
                -3.0 * gamma / (2.0 * (gamma + 1.0))
                + (s_grid / T) * 3.0 * gamma * (gamma + 2.0)
                / (2.0 * (gamma + 1.0) ** 2)
            )

            # Optional smoothing in s
            if smooth_sigma > 0.0:
                y = gaussian_filter1d(y, smooth_sigma)

            # Record minimum
            idx_min = np.argmin(y)
            s_min_list.append(s_grid[idx_min])
            y_min_list.append(y[idx_min])

        theory_xmins[T] = np.array(s_min_list)
        theory_ymins[T] = np.array(y_min_list)

    return theory_xmins, theory_ymins, s_grid, v0_values
