"""BENDS: Bayesian Estimation of Noisy Derivative Shifts."""

import numpy as np


def BENDS(t, y, t_N, t_D, PreL, PostL, tau, q):
    """
    Identify the candidate time near a documented event that produces the
    largest Bayesian-estimated change in signal slope in the expected
    direction.

    Parameters
    ----------
    t : array_like
        Relative time in minutes. Must be strictly increasing and uniformly
        sampled.
    y : array_like
        Signal to analyze. BENDS performs no preprocessing.
    t_N : float
        Documented event time in minutes.
    t_D : float
        Maximum timing uncertainty around t_N in minutes.
    PreL : float
        Pre-candidate window length in minutes.
    PostL : float
        Post-candidate window length in minutes.
    tau : float
        Positive residual-fluctuation tolerance.
    q : int
        Expected derivative-shift direction, either +1 or -1.

    Returns
    -------
    c_star, m1, m2, b1, b2 : tuple of floats
        Inferred onset time, pre/post Bayesian slopes, and pre/post line
        intercepts in the original signal scale.

    Notes
    -----
    The local coordinate for b1 starts at c_star - PreL.
    The local coordinate for b2 starts at c_star.
    """

    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)

    _validate_inputs(t, y, t_N, t_D, PreL, PostL, tau, q)

    # Candidate effect-onset times, C, evaluated at 1-min intervals.
    C = np.arange(t_N - t_D, t_N + t_D + 0.5, 1.0)

    best_Phi = -np.inf

    c_star = np.nan
    m1 = np.nan
    m2 = np.nan
    b1 = np.nan
    b2 = np.nan

    for c in C:

        # Window j = 1: pre-candidate window.
        idx1 = (t >= c - PreL) & (t <= c)
        x1 = t[idx1] - (c - PreL)
        y1 = y[idx1]

        # Window j = 2: post-candidate window.
        idx2 = (t >= c) & (t <= c + PostL)
        x2 = t[idx2] - c
        y2 = y[idx2]

        m1_c, b1_c = _bayesian_slope(x1, y1, PreL, tau)
        m2_c, b2_c = _bayesian_slope(x2, y2, PostL, tau)

        # Direction-adjusted derivative shift from Appendix B.
        Phi_c = q * (m2_c - m1_c)

        if Phi_c > best_Phi:
            best_Phi = Phi_c
            c_star = c
            m1 = m1_c
            m2 = m2_c
            b1 = b1_c
            b2 = b2_c

    return c_star, m1, m2, b1, b2


def _bayesian_slope(x_j, y_j, L_j, tau):
    """Bayesian slope estimator for one candidate window."""

    mu_j = np.mean(x_j)
    s_j = np.std(x_j, ddof=1)
    ybar_j = np.median(y_j)

    if ybar_j == 0:
        raise ValueError("A candidate window has zero median signal value.")

    z_j = (x_j - mu_j) / s_j
    w_j = y_j / ybar_j

    # Ordinary least-squares estimate and residuals.
    X_j = np.column_stack((np.ones(len(x_j)), z_j))
    beta_OLS_j = np.linalg.lstsq(X_j, w_j, rcond=None)[0]
    r_j = w_j - X_j @ beta_OLS_j

    # Effective residual variance from Appendix B.
    Delta_t = np.median(np.diff(x_j))
    sigma_eff_j_sq = (Delta_t / L_j) * np.sum(r_j**2)

    # MAP slope in normalized coordinates.
    beta1_MAP_j = beta_OLS_j[1] / (1.0 + sigma_eff_j_sq / tau**2)

    # z_j is centered, so the unpenalized MAP intercept is mean(w_j).
    beta0_MAP_j = np.mean(w_j)

    # Return slope and intercept to the original signal scale.
    m_j = (ybar_j / s_j) * beta1_MAP_j
    b_j = ybar_j * (beta0_MAP_j - beta1_MAP_j * mu_j / s_j)

    return float(m_j), float(b_j)


def _validate_inputs(t, y, t_N, t_D, PreL, PostL, tau, q):

    if t.size != y.size or t.size < 3:
        raise ValueError("t and y must have equal length and at least 3 samples.")

    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)):
        raise ValueError("t and y must contain finite values only.")

    dt = np.diff(t)

    if np.any(dt <= 0):
        raise ValueError("t must be strictly increasing.")

    Delta_t = np.median(dt)
    tol = max(1e-9, 1e-6 * abs(Delta_t))

    if np.any(np.abs(dt - Delta_t) > tol):
        raise ValueError("t must be uniformly sampled.")

    scalars = [t_N, t_D, PreL, PostL, tau]
    if not all(np.isscalar(v) and np.isfinite(v) for v in scalars):
        raise ValueError("t_N, t_D, PreL, PostL, and tau must be finite scalars.")

    if t_D < 0:
        raise ValueError("t_D must be nonnegative.")

    if PreL <= 0 or PostL <= 0:
        raise ValueError("PreL and PostL must be positive.")

    if tau <= 0:
        raise ValueError("tau must be positive.")

    if q not in (-1, 1):
        raise ValueError("q must be +1 or -1.")

    required_start = t_N - t_D - PreL
    required_end = t_N + t_D + PostL

    if t[0] > required_start + tol or t[-1] < required_end - tol:
        raise ValueError(
            "The input segment must cover the full interval "
            "[t_N - t_D - PreL, t_N + t_D + PostL]."
        )
