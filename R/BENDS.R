# BENDS: Bayesian Estimation of Noisy Derivative Shifts
#
# Identifies the candidate time near a documented event that produces the
# largest Bayesian-estimated change in signal slope in the expected direction.

BENDS <- function(t, y, t_N, t_D, PreL, PostL, tau, q) {

    t <- as.numeric(t)
    y <- as.numeric(y)

    .validate_inputs(t, y, t_N, t_D, PreL, PostL, tau, q)

    # Candidate effect-onset times, C, evaluated at 1-min intervals.
    C <- seq(from = t_N - t_D, to = t_N + t_D, by = 1)

    best_Phi <- -Inf

    c_star <- NA_real_
    m1 <- NA_real_
    m2 <- NA_real_
    b1 <- NA_real_
    b2 <- NA_real_

    for (c in C) {

        # Window j = 1: pre-candidate window.
        idx1 <- t >= c - PreL & t <= c
        x1 <- t[idx1] - (c - PreL)
        y1 <- y[idx1]

        # Window j = 2: post-candidate window.
        idx2 <- t >= c & t <= c + PostL
        x2 <- t[idx2] - c
        y2 <- y[idx2]

        fit1 <- .bayesian_slope(x1, y1, PreL, tau)
        fit2 <- .bayesian_slope(x2, y2, PostL, tau)

        # Direction-adjusted derivative shift from Appendix B.
        Phi_c <- q * (fit2$m - fit1$m)

        if (Phi_c > best_Phi) {
            best_Phi <- Phi_c
            c_star <- c
            m1 <- fit1$m
            m2 <- fit2$m
            b1 <- fit1$b
            b2 <- fit2$b
        }
    }

    return(list(
        c_star = c_star,
        m1 = m1,
        m2 = m2,
        b1 = b1,
        b2 = b2
    ))
}


.bayesian_slope <- function(x_j, y_j, L_j, tau) {
    # Bayesian slope estimator for one candidate window.
    # Variable names follow Appendix B of the manuscript.

    mu_j <- mean(x_j)
    s_j <- sd(x_j)
    ybar_j <- median(y_j)

    if (ybar_j == 0) {
        stop("A candidate window has zero median signal value.")
    }

    z_j <- (x_j - mu_j) / s_j
    w_j <- y_j / ybar_j

    # Ordinary least-squares estimate and residuals.
    X_j <- cbind(1, z_j)
    beta_OLS_j <- qr.solve(X_j, w_j)
    r_j <- as.numeric(w_j - X_j %*% beta_OLS_j)

    # Effective residual variance from Appendix B.
    Delta_t <- median(diff(x_j))
    sigma_eff_j_sq <- (Delta_t / L_j) * sum(r_j^2)

    # MAP slope in normalized coordinates.
    beta1_MAP_j <- beta_OLS_j[2] / (1 + sigma_eff_j_sq / tau^2)

    # z_j is centered, so the unpenalized MAP intercept is mean(w_j).
    beta0_MAP_j <- mean(w_j)

    # Return slope and intercept to the original signal scale.
    m_j <- (ybar_j / s_j) * beta1_MAP_j
    b_j <- ybar_j * (beta0_MAP_j - beta1_MAP_j * mu_j / s_j)

    return(list(
        m = as.numeric(m_j),
        b = as.numeric(b_j)
    ))
}


.validate_inputs <- function(t, y, t_N, t_D, PreL, PostL, tau, q) {

    if (length(t) != length(y) || length(t) < 3) {
        stop("t and y must have equal length and at least 3 samples.")
    }

    if (any(!is.finite(t)) || any(!is.finite(y))) {
        stop("t and y must contain finite values only.")
    }

    dt <- diff(t)

    if (any(dt <= 0)) {
        stop("t must be strictly increasing.")
    }

    Delta_t <- median(dt)
    tol <- max(1e-9, 1e-6 * abs(Delta_t))

    if (any(abs(dt - Delta_t) > tol)) {
        stop("t must be uniformly sampled.")
    }

    scalar_finite <- function(x) {
        length(x) == 1 && is.finite(x)
    }

    if (!scalar_finite(t_N) || !scalar_finite(t_D) ||
        !scalar_finite(PreL) || !scalar_finite(PostL) ||
        !scalar_finite(tau)) {
        stop("t_N, t_D, PreL, PostL, and tau must be finite scalars.")
    }

    if (t_D < 0) {
        stop("t_D must be nonnegative.")
    }

    if (PreL <= 0 || PostL <= 0) {
        stop("PreL and PostL must be positive.")
    }

    if (tau <= 0) {
        stop("tau must be positive.")
    }

    if (!(length(q) == 1 && q %in% c(-1, 1))) {
        stop("q must be +1 or -1.")
    }

    required_start <- t_N - t_D - PreL
    required_end <- t_N + t_D + PostL

    if (t[1] > required_start + tol ||
        t[length(t)] < required_end - tol) {
        stop(
            paste0(
                "The input segment must cover the full interval ",
                "[t_N - t_D - PreL, t_N + t_D + PostL]."
            )
        )
    }
}
