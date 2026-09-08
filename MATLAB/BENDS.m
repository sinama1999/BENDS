function [c_star, m1, m2, b1, b2] = BENDS(t, y, t_N, t_D, PreL, PostL, tau, q)
% BENDS Bayesian Estimation of Noisy Derivative Shifts.
%
% Identifies the candidate time near a documented event that produces the
% largest Bayesian-estimated change in signal slope in the expected
% direction.
%
% INPUTS
%   t       : Nx1 numeric time vector [min], strictly increasing and
%             uniformly sampled
%   y       : Nx1 numeric signal vector, already preprocessed as desired
%   t_N     : documented event time [min]
%   t_D     : maximum timing uncertainty around t_N [min]
%   PreL    : pre-candidate window length [min]
%   PostL   : post-candidate window length [min]
%   tau     : residual-fluctuation tolerance, positive scalar
%   q       : expected derivative-shift direction, +1 or -1
%
% OUTPUTS
%   c_star  : inferred effect-onset time [min]
%   m1      : Bayesian-estimated pre-candidate slope [signal units/min]
%   m2      : Bayesian-estimated post-candidate slope [signal units/min]
%   b1      : pre-candidate line intercept [signal units]
%   b2      : post-candidate line intercept [signal units]
%
% The local coordinate for b1 starts at c_star - PreL.
% The local coordinate for b2 starts at c_star.
%
% BENDS performs no preprocessing.

    t = t(:);
    y = y(:);

    validate_inputs(t, y, t_N, t_D, PreL, PostL, tau, q);

    % Candidate effect-onset times, C, evaluated at 1-min intervals.
    C = (t_N - t_D : 1 : t_N + t_D).';

    best_Phi = -inf;

    c_star = NaN;
    m1 = NaN;
    m2 = NaN;
    b1 = NaN;
    b2 = NaN;

    for k = 1:numel(C)

        c = C(k);

        % Window j = 1: pre-candidate window.
        idx1 = t >= c - PreL & t <= c;
        x1 = t(idx1) - (c - PreL);
        y1 = y(idx1);

        % Window j = 2: post-candidate window.
        idx2 = t >= c & t <= c + PostL;
        x2 = t(idx2) - c;
        y2 = y(idx2);

        [m1_c, b1_c] = bayesian_slope(x1, y1, PreL, tau);
        [m2_c, b2_c] = bayesian_slope(x2, y2, PostL, tau);

        % Direction-adjusted derivative shift from Appendix B.
        Phi_c = q * (m2_c - m1_c);

        if Phi_c > best_Phi
            best_Phi = Phi_c;
            c_star = c;
            m1 = m1_c;
            m2 = m2_c;
            b1 = b1_c;
            b2 = b2_c;
        end
    end
end


function [m_j, b_j] = bayesian_slope(x_j, y_j, L_j, tau)
% Bayesian slope estimator for one candidate window.
% Variable names follow Appendix B of the manuscript.

    mu_j = mean(x_j);
    s_j = std(x_j);
    ybar_j = median(y_j);

    if ybar_j == 0
        error('BENDS:ZeroMedian', ...
            'A candidate window has zero median signal value.');
    end

    z_j = (x_j - mu_j) / s_j;
    w_j = y_j / ybar_j;

    % Ordinary least-squares estimate and residuals.
    X_j = [ones(numel(x_j), 1), z_j];
    beta_OLS_j = X_j \ w_j;
    r_j = w_j - X_j * beta_OLS_j;

    % Effective residual variance from Appendix B.
    Delta_t = median(diff(x_j));
    sigma_eff_j_sq = (Delta_t / L_j) * sum(r_j.^2);

    % MAP slope in normalized coordinates.
    beta1_MAP_j = beta_OLS_j(2) / ...
        (1 + sigma_eff_j_sq / tau^2);

    % z_j is centered, so the unpenalized MAP intercept is mean(w_j).
    beta0_MAP_j = mean(w_j);

    % Return slope and intercept to the original signal scale.
    m_j = (ybar_j / s_j) * beta1_MAP_j;
    b_j = ybar_j * ...
        (beta0_MAP_j - beta1_MAP_j * mu_j / s_j);
end


function validate_inputs(t, y, t_N, t_D, PreL, PostL, tau, q)

    if ~isnumeric(t) || ~isnumeric(y) || numel(t) ~= numel(y) || numel(t) < 3
        error('BENDS:InvalidData', ...
            't and y must be numeric vectors of equal length with at least 3 samples.');
    end

    if any(~isfinite(t)) || any(~isfinite(y))
        error('BENDS:InvalidData', ...
            't and y must contain finite values only.');
    end

    dt = diff(t);

    if any(dt <= 0)
        error('BENDS:InvalidTime', ...
            't must be strictly increasing.');
    end

    Delta_t = median(dt);
    tol = max(1e-9, 1e-6 * abs(Delta_t));

    if any(abs(dt - Delta_t) > tol)
        error('BENDS:NonuniformTime', ...
            't must be uniformly sampled.');
    end

    if ~isscalar(t_N) || ~isfinite(t_N)
        error('BENDS:InvalidParameter', 't_N must be a finite scalar.');
    end

    if ~isscalar(t_D) || ~isfinite(t_D) || t_D < 0
        error('BENDS:InvalidParameter', 't_D must be a nonnegative scalar.');
    end

    if ~isscalar(PreL) || ~isfinite(PreL) || PreL <= 0 || ...
       ~isscalar(PostL) || ~isfinite(PostL) || PostL <= 0
        error('BENDS:InvalidParameter', ...
            'PreL and PostL must be positive scalars.');
    end

    if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
        error('BENDS:InvalidParameter', 'tau must be a positive scalar.');
    end

    if ~isscalar(q) || ~ismember(q, [-1, 1])
        error('BENDS:InvalidDirection', 'q must be +1 or -1.');
    end

    required_start = t_N - t_D - PreL;
    required_end = t_N + t_D + PostL;

    if t(1) > required_start + tol || t(end) < required_end - tol
        error('BENDS:InsufficientCoverage', ...
            ['The input segment must cover the full interval ', ...
             '[t_N - t_D - PreL, t_N + t_D + PostL].']);
    end
end
