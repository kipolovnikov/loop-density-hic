import numpy as np
import os
import pickle

# -------------------------------------------------------------------
# Helper: load (slope, mids) and extract (s_min, y_min) for one state
# -------------------------------------------------------------------

def load_minimum(prefix: str, state: str, data_dir: str, n_points: int = 20):
    """
    Load precomputed log-derivative for a given dataset/state and
    extract the position and depth of the minimum (s_min, y_min)
    over the first `n_points` bins.

    Parameters
    ----------
    prefix : str
        Prefix used in the pickle filenames, e.g. "abramo", "rao", "bonev".
    state : str
        State / cooler-name suffix used in the filenames.
    data_dir : str
        Directory containing full_logder_x_*.pickle and full_logder_y_*.pickle.
    n_points : int
        Number of initial points in the curve to search for the minimum.

    Returns
    -------
    x_min_bp : float
        Genomic separation at the minimum, in base pairs.
    y_min : float
        Value of the log-derivative at the minimum.
    """
    y_path = os.path.join(data_dir, f"full_logder_y_{prefix}_{state}.pickle")
    x_path = os.path.join(data_dir, f"full_logder_x_{prefix}_{state}.pickle")

    with open(y_path, "rb") as f:
        slope = pickle.load(f)
    with open(x_path, "rb") as f:
        mids = pickle.load(f)

    # restrict search to first n_points (short-scale regime)
    y_slice = slope[:n_points]
    idx_min = np.argmin(y_slice)
    y_min = y_slice[idx_min]

    # mids are in kb; convert to bp (x_min_bp)
    x_min_kb = mids[1 + idx_min]
    x_min_bp = x_min_kb * 1000.0

    return x_min_bp, y_min


# -------------------------------------------------------------------
# Compute minima for each dataset
# -------------------------------------------------------------------

def collect_minima(prefix: str, states, data_dir: str, n_points: int = 20):
    
    xs, ys = [], []
    for state in states:
        x_min, y_min = load_minimum(prefix, state, data_dir, n_points=n_points)
        xs.append(x_min)
        ys.append(y_min)
    return np.array(xs), np.array(ys)


# -------------------------------------------------------------------
# Load the slope for each dataset
# -------------------------------------------------------------------

def load_data(prefix: str, state: str, data_dir: str):

    y_path = os.path.join(data_dir, f"full_logder_y_{prefix}_{state}.pickle")
    x_path = os.path.join(data_dir, f"full_logder_x_{prefix}_{state}.pickle")

    with open(y_path, "rb") as f:
        slope = pickle.load(f)
    with open(x_path, "rb") as f:
        mids = pickle.load(f)
        
    return mids, slope