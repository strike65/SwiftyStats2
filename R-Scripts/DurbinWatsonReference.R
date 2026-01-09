#!/usr/bin/env Rscript
# Reference calculations for Durbin-Watson statistics and Imhof p-values.
# Uses CompQuadForm::imhof for the quadratic form CDF and optionally lmtest::dwtest
# to cross-check the DW statistic on a simple regression model.

`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

suppressPackageStartupMessages({
  have_cqf <- requireNamespace("CompQuadForm", quietly = TRUE)
  have_lmtest <- requireNamespace("lmtest", quietly = TRUE)
})

if (!have_cqf) {
  stop("CompQuadForm package is required. Install with: install.packages('CompQuadForm')")
}

# Durbin-Watson statistic: sum(diff(residuals)^2) / sum(residuals^2)
dw_stat <- function(residuals) {
  if (length(residuals) < 2) return(NA_real_)
  num <- sum(diff(residuals)^2)
  den <- sum(residuals^2)
  if (den == 0) return(NA_real_)
  num / den
}

# Lower-tail CDF for Q = sum(lambda_i * Z_i^2) using Imhof.
imhof_cdf <- function(q, lambda, ...) {
  res <- CompQuadForm::imhof(q = q, lambda = lambda, ...)
  prob <- res$value %||% res$prob %||% res$Qq
  list(prob = prob, ifault = res$ifault)
}

# Durbin-Watson lower-tail p-value: lambda_i = eigenvalues - dObs, P(Q <= 0)
dw_p_lower <- function(dObs, eigenvalues, ...) {
  lambda <- eigenvalues - dObs
  imhof_cdf(q = 0, lambda = lambda, ...)$prob
}

cat("### Durbin-Watson reference calculations ###\n")

# Basic residuals (matches Swift unit test)
res_basic <- c(1.0, 2.0, 3.0)
dw_basic <- dw_stat(res_basic)
cat(sprintf("DW(res_basic) = %.15f\n", dw_basic))

# Longer residual vector (matches Swift stress test)
res_long <- c(0.25, -1.1, 3.4, -0.6, 2.7, -2.9, 0.8, -0.4, 1.3, -1.7, 0.5)
dw_long <- dw_stat(res_long)
cat(sprintf("DW(res_long)  = %.15f\n", dw_long))

# Optional cross-check against lmtest::dwtest on a simple AR(1)-style regression
if (have_lmtest) {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 0.5 * x + rnorm(n, sd = 0.75)
  fit <- lm(y ~ x)
  dw_lm <- lmtest::dwtest(fit)
  cat(sprintf("lmtest::dwtest statistic = %.6f (p-value = %.6f)\n", dw_lm$statistic[["DW"]], dw_lm$p.value))
} else {
  cat("lmtest not installed; skipping lmtest::dwtest cross-check.\n")
}

if (have_lmtest) {
  cat("\n### dwtest on an arbitrary numeric time series ###\n")
  ts_args <- commandArgs(trailingOnly = TRUE)
  parse_series <- function(args) {
    if (length(args) == 0) return(NULL)
    parts <- unlist(strsplit(paste(args, collapse = " "), "[, ]+"))
    vals <- as.numeric(parts)
    vals[is.finite(vals)]
  }
  series_values <- parse_series(ts_args)
  if (!is.null(series_values) && length(series_values) >= 4) {
    y <- series_values[-1]
    xlag <- series_values[-length(series_values)]
    df <- data.frame(y = y, xlag = xlag)
    fit <- lm(y ~ xlag, data = df)
    dw_res <- lmtest::dwtest(fit)
    cat(sprintf("Input series length: %d\n", length(series_values)))
    cat(sprintf("dwtest statistic    : %.6f\n", dw_res$statistic[["DW"]]))
    cat(sprintf("dwtest p-value      : %.6f\n", dw_res$p.value))
  } else {
    cat("Pass a numeric series via command args (comma/space separated) to run dwtest on y_t ~ y_{t-1}.\n")
  }
} else {
  cat("\nInstall 'lmtest' to run dwtest on custom series.\n")
}

cat("\n### Imhof quadratic form checks ###\n")

# Single lambda should match scaled Chi-square CDF: P(lambda * X^2 <= x) = P(ChiSq(1) <= x / lambda)
lambda_single <- 1.5
x_single <- 2.0
imhof_single <- imhof_cdf(q = x_single, lambda = lambda_single)$prob
chisq_single <- pchisq(x_single / lambda_single, df = 1)
cat(sprintf("Imhof single-lambda CDF (x=%.3f, lambda=%.3f): %.8f\n", x_single, lambda_single, imhof_single))
cat(sprintf("Chi-square reference CDF                      : %.8f\n", chisq_single))

# Symmetric lambdas imply P(Q <= 0) ≈ 0.5
lambda_sym <- c(1.0, -1.0)
p_sym <- imhof_cdf(q = 0, lambda = lambda_sym)$prob
cat(sprintf("Symmetric lambda CDF at 0 (lambda=+1,-1): %.4f\n", p_sym))

# Durbin-Watson lower-tail p-value for eigenvalues = 1.0, dObs = 1.0 (degenerate lambda=0)
p_dw_degenerate <- dw_p_lower(dObs = 1.0, eigenvalues = 1.0)
cat(sprintf("DW lower-tail p-value (eigen=1, dObs=1) : %.4f\n", p_dw_degenerate))

# Durbin-Watson lower-tail p-value for eigenvalues matching Swift chi-square test
p_dw_chisq <- dw_p_lower(dObs = 0.0, eigenvalues = lambda_single)
cat(sprintf("DW lower-tail p-value (lambda=1.5, x=2) : %.8f\n", p_dw_chisq))
