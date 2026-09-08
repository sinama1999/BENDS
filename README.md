# BENDS

**BENDS: Bayesian Estimation of Noisy Derivative Shifts**

BENDS identifies the time of an expected change in signal slope near a documented event time. It was developed for the BP-guided effect-onset inference method described in the manuscript *Blood Pressure-Guided Inference of Vasoactive Drug Effect Onset for More Reliable Dose-Response Analysis*.

The MATLAB, Python, and R implementations in this repository use the same calculations and the same input definitions.

## What BENDS does

Given:

- a uniformly sampled signal `y`,
- its relative time vector `t` in minutes,
- a nurse-documented event time `t_N`,
- a maximum timing uncertainty `t_D`,
- pre-candidate and post-candidate window lengths `PreL` and `PostL`,
- a residual-fluctuation tolerance `tau`, and
- an expected derivative-shift direction `q`,

BENDS evaluates candidate effect-onset times at 1-minute intervals from

` t_N - t_D ` to ` t_N + t_D `.

For every candidate time, the signal slope is estimated before and after the candidate using the Bayesian linear regression described in Appendix B of the manuscript. The selected onset time is the candidate that maximizes

`Phi(c) = q * (m2 - m1)`

where `m1` and `m2` are the estimated pre-candidate and post-candidate slopes.

Use:

- `q = +1` when the expected change is toward a more positive slope, meaning `m2 - m1 > 0`.
- `q = -1` when the expected change is toward a more negative slope, meaning `m2 - m1 < 0`.

BENDS estimates the onset of the observable signal change. It does not estimate the true physical event time unless those two times coincide.

## Inputs

All three implementations use the same input order:

```text
BENDS(t, y, t_N, t_D, PreL, PostL, tau, q)
```

| Input | Meaning | Units / values |
| --- | --- | --- |
| `t` | Relative time vector | minutes |
| `y` | Signal to analyze | original signal units |
| `t_N` | Nurse-documented event time | minutes |
| `t_D` | Maximum timing uncertainty around `t_N` | minutes |
| `PreL` | Pre-candidate regression-window length | minutes |
| `PostL` | Post-candidate regression-window length | minutes |
| `tau` | Residual-fluctuation tolerance used by the Bayesian slope estimator | positive scalar |
| `q` | Expected direction of the derivative shift | `+1` or `-1` |

### Required signal properties

The BENDS function performs no preprocessing. The supplied `y` vector must already be in the form that should be analyzed.

For exact reproduction of the manuscript method:

- `t` must be numeric, strictly increasing, uniformly sampled, and expressed in minutes.
- `y` must have the same length as `t` and contain finite values only.
- Missing samples should be handled before calling BENDS.
- Signal smoothing, interpolation, artifact handling, and any other preprocessing must be performed outside BENDS.
- In the manuscript application, BENDS was applied to the already preprocessed BP signal.
- Because each local signal window is normalized by its median, its median must be nonzero.
- The input segment must cover the complete search range and both regression windows:

```text
[t_N - t_D - PreL,  t_N + t_D + PostL]
```

The manuscript used `t_D = 20 min`, `PreL = 15 min`, `PostL = 15 min`, `tau = 0.05`, and a 1-minute candidate spacing. The candidate spacing is fixed at 1 minute in these implementations to reproduce the method.

## Outputs

### MATLAB and Python

```text
c_star, m1, m2, b1, b2
```

### R

The function returns a named list containing the same five quantities.

| Output | Meaning |
| --- | --- |
| `c_star` | BENDS-inferred effect-onset time, in minutes |
| `m1` | Bayesian-estimated slope in the pre-candidate window |
| `m2` | Bayesian-estimated slope in the post-candidate window |
| `b1` | Intercept of the fitted pre-candidate line in original signal units |
| `b2` | Intercept of the fitted post-candidate line in original signal units |

The local time coordinate for `b1` begins at the start of the pre-candidate window. The local time coordinate for `b2` begins at `c_star`. Therefore, the fitted lines are

```text
pre:  y = b1 + m1 * x1,    x1 = t - (c_star - PreL)
post: y = b2 + m2 * x2,    x2 = t - c_star
```

## Bayesian slope estimation

For each window `j`, BENDS follows the Appendix B formulation.

The local time and signal are transformed as

```text
z_j = (x_j - mu_j) / s_j
w_j = y_j / ybar_j
```

where `mu_j` and `s_j` are the mean and sample standard deviation of local time and `ybar_j` is the median signal value in that window.

Ordinary least squares is first used to obtain `beta_OLS,j` and residuals `r_j`. The effective residual variance is

```text
sigma_eff,j^2 = (Delta_t / L_j) * sum(r_j^2)
```

where `L_j` is `PreL` for the first window and `PostL` for the second window.

With the zero-centered Gaussian slope prior used in the manuscript, the MAP slope in normalized coordinates is

```text
beta1_MAP,j = beta1_OLS,j / (1 + sigma_eff,j^2 / tau^2)
```

and the slope in the original signal scale is

```text
m_j = (ybar_j / s_j) * beta1_MAP,j
```

The candidate objective is

```text
Phi(c) = q * (m2 - m1)
```

and BENDS returns

```text
c_star = argmax Phi(c)
```

## Synthetic example

`example_data/synthetic_segment.csv` contains a made-up signal segment with:

- sampling interval: 0.5 min
- documented event time: `t_N = 35 min`
- timing uncertainty: `t_D = 20 min`
- pre-window length: `PreL = 15 min`
- post-window length: `PostL = 15 min`
- `tau = 0.05`
- expected derivative-shift direction: `q = +1`
- known synthetic slope-change time: `30 min`

The segment contains a negative pre-change slope, a positive post-change slope, and deterministic oscillatory noise. No preprocessing is performed by the example scripts.

Each language-specific example loads the same CSV file, runs BENDS, prints the five outputs, and plots the signal together with the documented time, known synthetic change time, inferred time, and fitted pre/post lines.

For the included synthetic segment, the inferred time is `c_star = 30 min`.

## Repository structure

```text
BENDS/
├── README.md
├── example_data/
│   └── synthetic_segment.csv
├── MATLAB/
│   ├── BENDS.m
│   └── example_BENDS.m
├── Python/
│   ├── bends.py
│   └── example_bends.py
└── R/
    ├── BENDS.R
    └── example_BENDS.R
```

## Running the examples

### MATLAB

From the repository root:

```matlab
run("MATLAB/example_BENDS.m")
```

### Python

Requirements:

```text
numpy
matplotlib
```

From the repository root:

```bash
python Python/example_bends.py
```

### R

The BENDS implementation itself uses base R only. From the repository root:

```bash
Rscript R/example_BENDS.R
```
