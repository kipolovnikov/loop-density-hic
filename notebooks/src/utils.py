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


def get_full_theory_curves(
    T_values=None,
    v0_values=None,
    s_min=1.0,
    s_max=1e4,
    s_step=0.01,
    smooth_sigma=0.0,
):
    """
    Compute the full theoretical log-derivative curves y(s) for a grid of
    loop periods T and effective fragment lengths v0.

    The theory used here is Eq. (3) in the main text:
        y(s) = d log P(s) / d log s
             = -3 γ / [2 (γ + 1)]
               + (s / T) * 3 γ (γ + 2) / [2 (γ + 1)^2]
        with γ(s) = 2 s / (3 v0)

    Parameters
    ----------
    T_values : array-like, optional
        Loop periods T (in kb) to evaluate. If None, defaults to
        np.arange(100, 301, 10), i.e. 100, 110, ..., 300 kb.
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
        Step size for genomic separation grid (in kb).
        Default: 0.01
    smooth_sigma : float, optional
        Standard deviation of the Gaussian kernel for optional smoothing
        in index space (not in kb). If 0.0, no smoothing is applied.
        Default: 0.0

    Returns
    -------
    s_grid : np.ndarray, shape (n_s,)
        Genomic separation grid (in kb).
    y_array : np.ndarray, shape (n_T, n_v0, n_s)
        Theoretical log-derivative curves y(s) for each (T, v0) pair.
        Indexing: y_array[i_T, j_v0, :] corresponds to T_values[i_T],
        v0_values[j_v0].
    T_values : np.ndarray, shape (n_T,)
        Array of loop periods actually used (in kb).
    v0_values : np.ndarray, shape (n_v0,)
        Array of effective fragment lengths actually used (in kb).
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
    n_s = s_grid.size
    n_T = T_values.size
    n_v0 = v0_values.size

    # Allocate output: shape (n_T, n_v0, n_s)
    y_array = np.empty((n_T, n_v0, n_s), dtype=float)

    for i_T, T in enumerate(T_values):
        for j_v0, v0 in enumerate(v0_values):
            # γ(s) = 2s / (3 v0)
            gamma = 2.0 * s_grid / (3.0 * v0)

            # y(s) = -3γ / [2(γ+1)] + (s/T) * 3γ(γ+2) / [2(γ+1)^2]
            y = (
                -3.0 * gamma / (2.0 * (gamma + 1.0))
                + (s_grid / T)
                * 3.0 * gamma * (gamma + 2.0)
                / (2.0 * (gamma + 1.0) ** 2)
            )

            if smooth_sigma > 0.0:
                y = gaussian_filter1d(y, smooth_sigma)

            y_array[i_T, j_v0, :] = y

    return s_grid, y_array, T_values, v0_values




def get_best_params(
    mids,
    slope,
    s_grid,
    y_array,
    T_values,
    v0_values,
    sfit_min=3.0,
    sfit_max=50.0,
    sfit_step=1.0,
    return_curves=False,
):
    """
    Find best-fit (T, v0) by comparing an experimental log-derivative (slope)
    to a grid of theoretical curves y_array(T, v0, s).

    Parameters
    ----------
    mids : array-like
        Genomic separation midpoints (in kb) corresponding to `slope`.
        Typically mids[1:] are used for the log-derivative.
    slope : array-like
        Experimental log-derivative values d log P(s) / d log s
        (same length as mids[1:]).
    s_grid : array-like, shape (n_s,)
        Genomic separation grid (in kb) used for the theoretical curves.
    y_array : array-like, shape (n_T, n_v0, n_s)
        Theoretical log-derivative curves y(s) from `get_full_theory_curves`.
        y_array[i_T, j_v0, :] corresponds to T_values[i_T], v0_values[j_v0].
    T_values : array-like, shape (n_T,)
        Loop periods T (in kb) for the theoretical grid.
    v0_values : array-like, shape (n_v0,)
        Effective fragment lengths v0 (in kb) for the theoretical grid.
    sfit_min : float, optional
        Minimum genomic separation (in kb) used for fitting.
        Default: 3.0
    sfit_max : float, optional
        Maximum genomic separation (in kb) used for fitting.
        Default: 50.0
    sfit_step : float, optional
        Step (in kb) for the fitting grid.
        Default: 1.0
    return_curves : bool, optional
        If True, also return the experimental and best-fit theoretical
        curves interpolated onto the fitting grid. Default: True.

    Returns
    -------
    T_opt : float
        Best-fit loop period (in kb).
    v0_opt : float
        Best-fit effective fragment length (in kb).
    s_fit : np.ndarray, shape (n_fit,)
        Fitting genomic separation grid (in kb).
    y_exp_fit : np.ndarray, shape (n_fit,)   (if return_curves=True)
        Experimental log-derivative interpolated onto s_fit.
    y_th_fit : np.ndarray, shape (n_fit,)    (if return_curves=True)
        Best-fit theoretical log-derivative interpolated onto s_fit.
    """

    mids = np.asarray(mids, dtype=float)
    slope = np.asarray(slope, dtype=float)
    s_grid = np.asarray(s_grid, dtype=float)
    T_values = np.asarray(T_values, dtype=float)
    v0_values = np.asarray(v0_values, dtype=float)

    # Fitting grid in s (kb)
    s_fit = np.arange(sfit_min, sfit_max, sfit_step, dtype=float)

    # Interpolate experimental slope onto s_fit
    # (mids[1:] because slope is typically computed on the midpoints between bins)
    y_exp_fit = np.interp(s_fit, mids[1:], slope)

    best_err = np.inf
    T_opt = None
    v0_opt = None
    y_th_fit_opt = None

    # Loop over all (T, v0) combinations
    n_T = T_values.size
    n_v0 = v0_values.size

    for i_T, T in enumerate(T_values):
        for j_v0, v0 in enumerate(v0_values):

            # Interpolate theoretical curve onto s_fit
            y_th = np.interp(s_fit, s_grid, y_array[i_T, j_v0, :])

            # Simple least-squares error over the fitting window
            err = np.sum((y_th - y_exp_fit) ** 2)

            if err < best_err:
                best_err = err
                T_opt = T
                v0_opt = v0
                y_th_fit_opt = y_th

    if return_curves:
        return T_opt, v0_opt, s_fit, y_th_fit_opt
    else:
        return T_opt, v0_opt

