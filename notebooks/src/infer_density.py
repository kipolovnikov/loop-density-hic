#!/usr/bin/env python3

"""
infer_density.py

Infer loop period T and effective fragment length v0 from a single
Hi-C contact curve (either raw scaling P(s) or its log-derivative),
and plot the experimental log-derivative together with the best-fit
theoretical curve.

Usage examples
--------------
# If y-file already contains log-derivative (slope)
python src/infer_density.py \
    --x data/full_logder_x_rao_GM12878.pickle \
    --y data/full_logder_y_rao_GM12878.pickle \
    --mode slope \
    --output-plot results/gm12878_fit.png

# If y-file contains raw scaling P(s)
python src/infer_density.py \
    --x data/scaling_mids.pickle \
    --y data/scaling_values.pickle \
    --mode scaling
"""

import argparse
import os
import pickle

import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Theory utilities
# ----------------------------------------------------------------------

def theory_logderivative(s_kb, T_kb, v0_kb, smooth_sigma=0.1):
    """
    Logarithmic derivative y(s) for the looped polymer with finite capture radius.

    Parameters
    ----------
    s_kb : array-like
        Genomic separation in kilobases.
    T_kb : float
        Loop period (kb).
    v0_kb : float
        Effective fragment length (kb).
    smooth_sigma : float
        Gaussian smoothing sigma (in index units) applied to the curve.

    Returns
    -------
    y : np.ndarray
        Theoretical log-derivative y(s).
    """
    s_kb = np.asarray(s_kb, dtype=float)
    gamma = 2.0 * s_kb / (3.0 * v0_kb)

    # Eq. (3) from the main text:
    # y(s) = -3γ / [2(γ+1)] + (s/T) * 3γ(γ+2) / [2(γ+1)^2]
    first_term = -3.0 * gamma / (2.0 * (gamma + 1.0))
    second_term = (s_kb / T_kb) * 3.0 * gamma * (gamma + 2.0) / (2.0 * (gamma + 1.0) ** 2)
    y = first_term + second_term

    if smooth_sigma is not None and smooth_sigma > 0:
        y = gaussian_filter1d(y, smooth_sigma)

    return y


def compute_theory_minima(T_values, v0_values,
                          s_min_kb=1.0, s_max_kb=1000.0, ds_kb=0.01,
                          smooth_sigma=0.1):
    """
    Precompute (s_min, y_min) for a grid of (T, v0) theory curves.

    Parameters
    ----------
    T_values : sequence of float
        Loop periods (kb) to scan.
    v0_values : sequence of float
        Effective fragment lengths (kb) to scan.
    s_min_kb, s_max_kb : float
        Range of s (kb) over which to compute the curves.
    ds_kb : float
        Step in s (kb).
    smooth_sigma : float
        Gaussian smoothing for the theoretical curves.

    Returns
    -------
    minima : dict
        Dictionary keyed by (T, v0) with values (s_min, y_min),
        both in kb and log-derivative units respectively.
    """
    s = np.arange(s_min_kb, s_max_kb, ds_kb)
    minima = {}

    for T in T_values:
        for v0 in v0_values:
            y = theory_logderivative(s, T_kb=T, v0_kb=v0, smooth_sigma=smooth_sigma)
            idx_min = np.argmin(y)
            s_min = s[idx_min]
            y_min = y[idx_min]
            minima[(T, v0)] = (s_min, y_min)

    return minima


def fit_T_v0_from_minimum(s_exp_kb, y_exp,
                          minima_dict,
                          x_weight=1e-4):
    """
    Fit T and v0 using only the position and depth of the minimum.

    We compare the experimental (s_min, y_min) to precomputed theory
    minima on a grid and select the (T, v0) pair with minimal
    weighted squared distance:
        err = x_weight * (s_th - s_exp)^2 + (y_th - y_exp)^2

    Parameters
    ----------
    s_exp_kb : float
        Experimental minimum position (kb).
    y_exp : float
        Experimental minimum value of the log-derivative.
    minima_dict : dict
        Output of compute_theory_minima.
    x_weight : float
        Weight given to differences in s_min relative to y_min.

    Returns
    -------
    T_best, v0_best, err_best
    """
    best_key = None
    best_err = np.inf

    for (T, v0), (s_th, y_th) in minima_dict.items():
        err = x_weight * (s_th - s_exp_kb) ** 2 + (y_th - y_exp) ** 2
        if err < best_err:
            best_err = err
            best_key = (T, v0)

    if best_key is None:
        raise RuntimeError("No minimum found in the theory grid (minima_dict is empty?).")

    T_best, v0_best = best_key
    return T_best, v0_best, best_err


# ----------------------------------------------------------------------
# Data utilities
# ----------------------------------------------------------------------

def load_array(path):
    """
    Load a numpy array or list from a pickle or .npy file.

    Parameters
    ----------
    path : str
        Path to file.

    Returns
    -------
    arr : np.ndarray
    """
    ext = os.path.splitext(path)[1]
    if ext == ".npy":
        return np.load(path)
    else:
        # default: pickle
        with open(path, "rb") as f:
            arr = pickle.load(f)
        return np.asarray(arr)


def compute_log_derivative_from_scaling(mids_kb, scaling, smooth_sigma=1.0):
    """
    Convert a raw scaling curve P(s) into a log-derivative.

    slope(s) = d log P / d log s, computed as discrete differences
    and optionally smoothed with a Gaussian filter.

    Parameters
    ----------
    mids_kb : array-like
        Genomic separations (kb).
    scaling : array-like
        Contact probability P(s).
    smooth_sigma : float
        Sigma for gaussian_filter1d applied to the slope.

    Returns
    -------
    slope : np.ndarray
        Log-derivative evaluated at mids_kb[1:].
    """
    mids_kb = np.asarray(mids_kb, dtype=float)
    scaling = np.asarray(scaling, dtype=float)

    # discrete derivative in log-log space
    slope = np.diff(np.log(scaling)) / np.diff(np.log(mids_kb))

    if smooth_sigma is not None and smooth_sigma > 0:
        slope = gaussian_filter1d(slope, smooth_sigma)

    return slope


def extract_minimum_from_slope(mids_kb, slope, n_points=20):
    """
    Extract (s_min, y_min) from the first n_points of the log-derivative.

    Follows the convention used in the paper:
    slope[i] is plotted at mids[i+1], so the minimum at slope[idx]
    corresponds to mids[idx+1].

    Parameters
    ----------
    mids_kb : array-like
        Genomic separations (kb).
    slope : array-like
        Log-derivative values.
    n_points : int
        Number of initial indices in slope to search over.

    Returns
    -------
    s_min_kb : float
        Position of the minimum in kb.
    y_min : float
        Minimum value of the log-derivative.
    """
    slope = np.asarray(slope, dtype=float)
    mids_kb = np.asarray(mids_kb, dtype=float)

    end = min(n_points, len(slope))
    slice_ = slope[:end]

    idx_min = np.argmin(slice_)
    y_min = slice_[idx_min]
    s_min_kb = mids_kb[idx_min + 1]  # +1 because slope[i] corresponds to mids[i+1]

    return s_min_kb, y_min


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def plot_fit(mids_kb, slope_exp, T_best, v0_best, output_path=None, show=True):
    """
    Plot experimental log-derivative and best-fit theoretical curve.

    Parameters
    ----------
    mids_kb : array-like
        Genomic separations (kb).
    slope_exp : array-like
        Experimental log-derivative values (same length as mids_kb[1:]).
    T_best : float
        Best-fit loop period (kb).
    v0_best : float
        Best-fit effective fragment length (kb).
    output_path : str or None
        If not None, save the figure to this path.
    show : bool
        If True, display the figure interactively.
    """
    mids_kb = np.asarray(mids_kb, dtype=float)
    slope_exp = np.asarray(slope_exp, dtype=float)

    # Theory evaluated at the same s as experimental slope
    s_for_theory = mids_kb[1:]
    slope_th = theory_logderivative(s_for_theory, T_kb=T_best, v0_kb=v0_best, smooth_sigma=0.1)

    fig = plt.figure(figsize=(5, 5)) 

    plt.plot(s_for_theory, slope_exp, color="k", lw=2, label="Experiment")
    plt.plot(s_for_theory, slope_th, color="r", lw=2, ls="--",
            label=f"Theory: T={T_best:.0f} kb, v0={v0_best:.2f} kb")

    plt.xlabel("Genomic separation s (kb)")
    plt.ylabel("log-derivative of P(s)")
    plt.legend()

    plt.ylim([-1.5, -0.2])
    plt.xlim([5, 100])

    if output_path is not None:
        fig.savefig(output_path, dpi=300)

    if show:
        plt.show()
    else:
        plt.close(fig)

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
    
# ----------------------------------------------------------------------
# Main CLI
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Infer loop density (T) and effective fragment length (v0) "
                    "from a single Hi-C curve and plot the best-fit theory."
    )
    parser.add_argument("--x", required=True,
                        help="Path to mids array (kb) [pickle or .npy].")
    parser.add_argument("--y", required=True,
                        help="Path to y-array: raw scaling P(s) or log-derivative.")
    parser.add_argument("--mode", choices=["scaling", "slope"], default="slope",
                        help="'scaling' if y is P(s); 'slope' if y is already log-derivative.")
    parser.add_argument("--smooth-sigma", type=float, default=1.0,
                        help="Gaussian sigma for scaling->slope conversion (default: 1.0).")
    parser.add_argument("--n-min-points", type=int, default=20,
                        help="Number of initial points in slope used to locate the minimum.")
    parser.add_argument("--sfit-min", type=int, default=8,
                        help="The left end of the genomic interval (in kb) where to fit the slope")
    parser.add_argument("--sfit-max", type=int, default=40,
                        help="The right end of the genomic interval (in kb) where to fit the slope")
    parser.add_argument("--output-plot", type=str, default=None,
                        help="Optional path to save the fit plot (PNG/PDF/etc.).")
    parser.add_argument("--no-show", action="store_true",
                        help="Do not display the plot interactively.")

    args = parser.parse_args()

    # Load data
    mids_kb = load_array(args.x)
    y_arr = load_array(args.y)

    # Compute log-derivative if necessary
    if args.mode == "scaling":
        slope = compute_log_derivative_from_scaling(mids_kb, y_arr,
                                                    smooth_sigma=args.smooth_sigma)
    else:
        slope = y_arr  # already a log-derivative
        # sanity check: expected length len(mids) - 1
        if len(slope) == len(mids_kb):
            print("[warning] slope has same length as mids; "
                  "usually it should be len(mids)-1.")

    # Extract experimental minimum
    s_min_kb, y_min = extract_minimum_from_slope(mids_kb, slope,
                                                 n_points=args.n_min_points)
    print(f"Experimental minimum: s_min = {s_min_kb:.2f} kb, y_min = {y_min:.3f}")

    # Precompute theory minima on a grid
    T_values = np.arange(100.0, 301.0, 10.0)  # 100–300 kb in 10 kb steps
    v0_values = np.arange(0.1, 15, 0.1)
    
    print("Precomputing theory minima grid ...")
    minima_dict = compute_theory_minima(
        T_values=T_values,
        v0_values=v0_values,
        s_min_kb=1.0,
        s_max_kb=1000.0,
        ds_kb=0.01,
        smooth_sigma=0.1,
    )

    # Fit T and v0
    T_best, v0_best, err_best = fit_T_v0_from_minimum(
        s_exp_kb=s_min_kb,
        y_exp=y_min,
        minima_dict=minima_dict,
        x_weight=1e-4,
    )

    print(f"Best-fit parameters from the DIP:")
    print(f"  T   = {T_best:.1f} kb  (loop period)")
    print(f"  v0  = {v0_best:.2f} kb (effective fragment length)")
    print(f"  density = {1000/T_best:.2f} loops/Mb")

    # Obtain the best-fit parameters from the CURVE around the dip, i.e. s \in (sfit_min, sfit_max)
    s_grid, y_array, T_values, v0_values = get_full_theory_curves(T_values, v0_values, s_min=1, s_max=50)

    T_opt, v0_opt = get_best_params(mids_kb, slope, s_grid, y_array, \
                        T_values, v0_values, sfit_min=args.sfit_min, sfit_max=args.sfit_max, return_curves=False) 

    print(f"Best fit from the CURVE: T = {T_opt:.1f} kb, v0 = {v0_opt:.2f} kb, density = {1000/T_opt:.2f} loops/Mb")


    # Plot experimental slope and best-fit theory
    plot_fit(
        mids_kb=mids_kb,
        slope_exp=slope,
        T_best=T_opt,
        v0_best=v0_opt,
        output_path=args.output_plot,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()
