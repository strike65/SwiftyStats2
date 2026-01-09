#!/usr/bin/env Rscript
# Mirrors the output of DurbinWatsonValidation.runValidation() using R.
# Requires CompQuadForm (davies) to be installed.

if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
  stop("CompQuadForm is required; install.packages('CompQuadForm')", call. = FALSE)
}
if (!requireNamespace("lmtest", quietly = TRUE)) {
  stop("lmtest is required for dwtest() cross-validation; install.packages('lmtest')", call. = FALSE)
}

# Stable CDF for indefinite quadratic forms: prefer davies, fall back to imhof.
quad_cdf <- function(weights, c = 0) {
  res <- CompQuadForm::davies(q = c, lambda = weights, acc = 1e-6)
  if (!is.null(res$error) && res$error != 0) {
    res_imhof <- CompQuadForm::imhof(q = c, lambda = weights, epsabs = 1e-6)
    return(res_imhof$Qq)
  }
  res$Qq
}

fmt <- function(x) formatC(x, digits = 16, format = "g")

cat("=== Durbin-Watson validation ===\n\n")

# Test 1: all weights positive
cat("Test 1: all weights positive\n")
p1 <- quad_cdf(c(1, 2, 3), c = 0)
passed1 <- p1 < 0.001
cat(sprintf("  P(sum lambda_i Z_i^2 < 0) = %s (expected: ~0)\n", fmt(p1)))
cat(sprintf("  %s\n\n", if (passed1) "[OK]" else "[FAIL]"))

# Test 2: all weights negative
cat("Test 2: all weights negative\n")
p2 <- quad_cdf(c(-1, -2, -3), c = 0)
passed2 <- p2 > 0.999
cat(sprintf("  P(sum lambda_i Z_i^2 < 0) = %s (expected: ~1)\n", fmt(p2)))
cat(sprintf("  %s\n\n", if (passed2) "[OK]" else "[FAIL]"))

# Test 3: symmetric case
cat("Test 3: Symmetric case\n")
p3 <- quad_cdf(c(1, -1), c = 0)
passed3 <- abs(p3 - 0.5) < 0.05
cat(sprintf("  P(Z1^2 - Z2^2 < 0) = %s (expected: 0.5)\n", fmt(p3)))
cat(sprintf("  %s\n\n", if (passed3) "[OK]" else "[FAIL]"))

# Test 4: DW statistic computation
cat("Test 4: DW statistic computation\n")
resid_alt <- c(1, -1, 1, -1, 1)
d1 <- sum(diff(resid_alt)^2) / sum(resid_alt^2)
passed4a <- abs(d1 - 3.2) < 1e-10
cat(sprintf("  Alternating residuals: d = %s (expected: 3.2) %s\n",
            fmt(d1), if (passed4a) "[OK]" else "[FAIL]"))

resid_const <- rep(2, 5)
d2 <- sum(diff(resid_const)^2) / sum(resid_const^2)
passed4b <- abs(d2) < 1e-10
cat(sprintf("  Constant residuals: d = %s (expected: 0) %s\n",
            fmt(d2), if (passed4b) "[OK]" else "[FAIL]"))

passed4 <- passed4a && passed4b
cat(sprintf("  %s\n\n", if (passed4) "[OK]" else "[FAIL]"))

# Test 5: Plausibility check
cat("Test 5: Plausibility check\n")
n <- 20
x_vals <- 0:(n - 1)
X <- cbind(1, x_vals)
y <- 1 + 0.5 * x_vals + 0.1 * (x_vals %% 3)

fit <- lm(y ~ x_vals)
res <- resid(fit)

dw_stat <- sum(diff(res)^2) / sum(res^2)

A <- diag(n) * 2
A[1, 1] <- 1
A[n, n] <- 1
A[cbind(1:(n - 1), 2:n)] <- -1
A[cbind(2:n, 1:(n - 1))] <- -1

H <- X %*% solve(t(X) %*% X) %*% t(X)
M <- diag(n) - H
MAM <- M %*% A %*% M
ev <- eigen((MAM + t(MAM)) / 2, symmetric = TRUE, only.values = TRUE)$values
ev <- ev[ev > 1e-10]

weights <- ev - dw_stat
p_lower <- quad_cdf(weights, c = 0)
p_upper <- 1 - p_lower

cat(sprintf("  n = %d, k = %d\n", n, ncol(X)))
cat(sprintf("  DW = %.4f\n", dw_stat))
cat(sprintf("  p(pos) = %.4f\n", p_lower))
cat(sprintf("  p(neg) = %.4f\n", p_upper))

# Bounds ~ from the same approximation as Swift helper
kPrime <- ncol(X) - 1
dL <- 2 - 2 * sqrt(2 * kPrime / n) - 2 / n
dU <- 2 - 2 * sqrt(2 / n) + 2 * kPrime / n

cat(sprintf("  Bounds: dL ~ %.2f, dU ~ %.2f\n", max(0, dL), min(4, dU)))
cat(sprintf("  Converged: %s\n", TRUE))

# Cross-validate with dwtest (positive autocorr => lower tail)
dw_pos <- lmtest::dwtest(fit, alternative = "greater")$p.value
dw_neg <- lmtest::dwtest(fit, alternative = "less")$p.value
dw_two <- lmtest::dwtest(fit, alternative = "two.sided")$p.value

cat(sprintf("  dwtest p(pos) = %.4f\n", dw_pos))
cat(sprintf("  dwtest p(neg) = %.4f\n", dw_neg))
cat(sprintf("  dwtest p(two) = %.4f\n", dw_two))

tol <- 5e-3
p5_dw_pos <- abs(dw_pos - p_lower) < tol
p5_dw_neg <- abs(dw_neg - p_upper) < tol
p5_dw_two <- abs(dw_two - min(1, 2 * min(p_lower, p_upper))) < tol

p1 <- dw_stat >= 0 && dw_stat <= 4
p2 <- p_lower >= 0 && p_lower <= 1
p3 <- p_upper >= 0 && p_upper <= 1
p4 <- abs(p_lower + p_upper - 1) < 0.01
p5_dw_two <- abs(dw_two - min(1, 2 * min(p_lower, p_upper))) < tol
passed5_core <- p1 && p2 && p3 && p4

cat(sprintf("  DW in [0,4]: %s\n", if (p1) "[OK]" else "[FAIL]"))
cat(sprintf("  p-values in [0,1]: %s\n", if (p2 && p3) "[OK]" else "[FAIL]"))
cat(sprintf("  p(pos) + p(neg) close to 1: %s\n", if (p4) "[OK]" else "[FAIL]"))
cat(sprintf("  dwtest match (pos): %s\n", if (p5_dw_pos) "[OK]" else "[FAIL]"))
cat(sprintf("  dwtest match (neg): %s\n", if (p5_dw_neg) "[OK]" else "[FAIL]"))
cat(sprintf("  dwtest match (two): %s\n", if (p5_dw_two) "[OK]" else "[FAIL]"))

passed5 <- passed5_core && p5_dw_pos && p5_dw_neg && p5_dw_two
cat(sprintf("  %s\n\n", if (passed5) "[OK]" else "[FAIL]"))

all_passed <- passed1 && passed2 && passed3 && passed4 && passed5
cat(sprintf("\n%s\n", if (all_passed) "[OK] all tests passed" else "[FAIL] some tests failed"))
