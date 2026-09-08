# BENDS synthetic example
# This example uses the same made-up segment as the MATLAB and Python examples.
# The input signal is used directly. BENDS performs no preprocessing.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
    script_dir <- dirname(script_path)
    repo_root <- dirname(script_dir)
} else {
    repo_root <- normalizePath(getwd())
    if (!file.exists(file.path(repo_root, "example_data", "synthetic_segment.csv"))) {
        repo_root <- dirname(repo_root)
    }
}

source(file.path(repo_root, "R", "BENDS.R"))
D <- read.csv(file.path(repo_root, "example_data", "synthetic_segment.csv"))

t <- D$t_min
y <- D$signal

# BENDS settings for this example.
t_N <- 35       # documented event time [min]
t_D <- 20       # maximum timing uncertainty [min]
PreL <- 15      # pre-candidate window length [min]
PostL <- 15     # post-candidate window length [min]
tau <- 0.05     # residual-fluctuation tolerance
q <- +1         # expected shift toward a more positive slope

# Known slope-change time used only to evaluate the synthetic example.
c_true <- 30

result <- BENDS(t, y, t_N, t_D, PreL, PostL, tau, q)

cat(sprintf("c_star = %.6f min\n", result$c_star))
cat(sprintf("m1     = %.6f signal units/min\n", result$m1))
cat(sprintf("m2     = %.6f signal units/min\n", result$m2))
cat(sprintf("b1     = %.6f signal units\n", result$b1))
cat(sprintf("b2     = %.6f signal units\n", result$b2))

# Fitted lines at the selected candidate.
idx_pre <- t >= result$c_star - PreL & t <= result$c_star
t_pre <- t[idx_pre]
x_pre <- t_pre - (result$c_star - PreL)
y_pre_fit <- result$b1 + result$m1 * x_pre

idx_post <- t >= result$c_star & t <= result$c_star + PostL
t_post <- t[idx_post]
x_post <- t_post - result$c_star
y_post_fit <- result$b2 + result$m2 * x_post

plot(
    t, y,
    type = "l",
    lwd = 1.2,
    xlab = "Time [min]",
    ylab = "Signal",
    main = "BENDS synthetic example"
)
lines(t_pre, y_pre_fit, lwd = 2)
lines(t_post, y_post_fit, lwd = 2)
abline(v = t_N, lty = 2)
abline(v = c_true, lty = 3)
abline(v = result$c_star, lty = 4)
grid()
legend(
    "topleft",
    legend = c(
        "Signal",
        "Pre-candidate fit",
        "Post-candidate fit",
        "Documented time",
        "Known synthetic change",
        "BENDS inferred time"
    ),
    lty = c(1, 1, 1, 2, 3, 4),
    lwd = c(1.2, 2, 2, 1, 1, 1),
    bty = "n"
)
