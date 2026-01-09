## KL divergence D_KL(Arcsine || N(0,1))
## Arcsine on [0,1]: p(x) = 1 / (pi * sqrt(x * (1 - x)))
## Standard normal:  q(x) = (1/sqrt(2*pi)) * exp(-x^2/2)

arcsine_pdf <- function(x) {
  inside <- (x > 0) & (x < 1)
  out <- numeric(length(x))
  out[inside] <- 1 / (pi * sqrt(x[inside] * (1 - x[inside])))
  out
}

normal_pdf <- function(x) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * x^2)
}

## Analytic closed form:
## D_KL = (2.5 * log(2)) - 0.5 * log(pi) + 3/16
kl_arcsine_to_stdnormal_analytic <- function() {
  (2.5 * log(2)) - 0.5 * log(pi) + 3/16
}

## Numerical check via integrate on [0,1]
kl_arcsine_to_stdnormal_numeric <- function() {
  integrand <- function(x) {
    p <- arcsine_pdf(x)
    q <- normal_pdf(x)
    p * (log(p) - log(q))
  }
  val <- integrate(integrand, lower = 0, upper = 1,
                   subdivisions = 1e6, rel.tol = 1e-10, abs.tol = 0)$value
  val
}

## Run
kl_analytic <- kl_arcsine_to_stdnormal_analytic()
kl_numeric  <- kl_arcsine_to_stdnormal_numeric()

cat(sprintf("KL(Arcsine || N(0,1)) analytic : %.12f\n", kl_analytic))
cat(sprintf("KL(Arcsine || N(0,1)) numeric  : %.12f\n", kl_numeric))

###
## KL divergence for continuous distributions in R (base only)
## We compute D_KL(P || Q) = ∫ p(x) [log p(x) - log q(x)] dx
## For Chi-square(k=22) vs Normal(mu=22, sd=sqrt(44)) on [0, +Inf).

## Generic KL using log-densities to improve numerical stability
kl_continuous <- function(logp, logq, lower, upper, ...) {
  integrand <- function(x) {
    lp <- logp(x, ...)
    lq <- logq(x, ...)
    p  <- exp(lp)
    p * (lp - lq)
  }
  integrate(integrand, lower = lower, upper = upper,
            rel.tol = 1e-9, abs.tol = 0)$value
}

## Log-pdf of Chi-square(k=22)
logpdf_chisq22 <- function(x) dchisq(x, df = 22, log = TRUE)

## Log-pdf of Normal with matching mean/variance: mu = 22, sd = sqrt(44)
logpdf_norm_match <- function(x) dnorm(x, mean = 22, sd = sqrt(44), log = TRUE)

## Optional: log-pdf of standard normal N(0,1)
logpdf_norm_std <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)

## Compute D_KL(ChiSq22 || Normal(22, sqrt(44)))
kl_chisq_vs_norm_match <- kl_continuous(
  logp  = logpdf_chisq22,
  logq  = logpdf_norm_match,
  lower = 0, upper = Inf
)

## (Optional) Compute D_KL(ChiSq22 || N(0,1)); finite but typically large
kl_chisq_vs_norm_std <- kl_continuous(
  logp  = logpdf_chisq22,
  logq  = logpdf_norm_std,
  lower = 0, upper = Inf
)

## Important note:
## D_KL(Normal || Chi-square) is infinite because Chi-square has support [0, +Inf)
## while Normal has positive mass on (-Inf, 0). Hence we do NOT compute Q||P.

cat(sprintf("KL(ChiSq(22) || N(22, sqrt(44))) = %.12f\n", kl_chisq_vs_norm_match))
cat(sprintf("KL(ChiSq(22) || N(0,1))         = %.12f\n", kl_chisq_vs_norm_std))

## KL divergence: Normal(mu, sd) || Cauchy(loc=2.5, scale=2.3)
## No extra packages required.

## Generic KL for continuous distributions using log-densities
kl_continuous <- function(logp, logq, lower = -Inf, upper = Inf, ...) {
  integrand <- function(x) {
    lp <- logp(x, ...)
    lq <- logq(x, ...)
    exp(lp) * (lp - lq)
  }
  integrate(integrand, lower = lower, upper = upper,
            rel.tol = 1e-9, abs.tol = 0)$value
}

## Log-pdf of Normal(mu, sd)
logpdf_normal <- function(x, mu, sd) dnorm(x, mean = mu, sd = sd, log = TRUE)

## Log-pdf of Cauchy(loc, scale)
logpdf_cauchy <- function(x, loc, scale) dcauchy(x, location = loc, scale = scale, log = TRUE)

## Wrapper: KL( Normal(mu, sd) || Cauchy(loc=2.5, scale=2.3) )
kl_normal_vs_cauchy <- function(mu = 0, sd = 1, loc = 2.5, scale = 2.3) {
  kl_continuous(
    logp = function(x) logpdf_normal(x, mu = mu, sd = sd),
    logq = function(x) logpdf_cauchy(x, loc = loc, scale = scale),
    lower = -Inf, upper = Inf
  )
}

## Optional: Monte-Carlo cross-check (importance sampling under P=N)
kl_mc_normal_vs_cauchy <- function(mu = 0, sd = 1, loc = 2.5, scale = 2.3, n = 1e6, seed = 123) {
  set.seed(seed)
  x <- rnorm(n, mean = mu, sd = sd)
  lp <- dnorm(x, mean = mu, sd = sd, log = TRUE)
  lq <- dcauchy(x, location = loc, scale = scale, log = TRUE)
  mean(lp - lq)
}

## Run: default compares N(0,1) to Cauchy(2.5, 2.3)
mu <- 0; sd <- 1
loc <- 2.5; scale <- 2.3

kl_int <- kl_normal_vs_cauchy(mu = mu, sd = sd, loc = loc, scale = scale)
kl_est <- kl_mc_normal_vs_cauchy(mu = mu, sd = sd, loc = loc, scale = scale, n = 2e5)

cat(sprintf("KL(N(%g,%g) || Cauchy(%g,%g)) [integral] : %.12f\n", mu, sd, loc, scale, kl_int))
cat(sprintf("KL(N(%g,%g) || Cauchy(%g,%g)) [MC est.]  : %.12f\n", mu, sd, loc, scale, kl_est))

## Note: D_KL(Cauchy || Normal) diverges to +Inf due to heavy tails of Cauchy.
