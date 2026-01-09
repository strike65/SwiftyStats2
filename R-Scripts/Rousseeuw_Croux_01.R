## ============================================================
## Monte Carlo estimation of finite-sample correction factors
## for Qn and Sn under N(0,1)
## Results are saved into the same directory as this script:
##   qn_dn_factors.txt : "n dn_hat" per line
##   sn_cn_factors.txt : "n cn_hat" per line
## ============================================================

set.seed(12345)

## ---------- Helper: determine script directory ----------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  if (length(idx) > 0) {
    script_path <- sub(file_arg, "", args[idx[1]])
    return(dirname(normalizePath(script_path)))
  } else {
    ## Fallback: current working directory (e.g. interactive run)
    return(getwd())
  }
}

script_dir <- get_script_dir()
cat("Script directory:", script_dir, "\n")

## ---------- Definition-based Qn and Sn (raw, no finite corr) ----------

Qn_raw <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (n < 2L) return(0)
  
  # pairwise distances |x_i - x_j|, i<j
  d <- numeric(n * (n - 1L) / 2L)
  k <- 1L
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      d[k] <- abs(x[i] - x[j])
      k <- k + 1L
    }
  }
  d <- sort(d)
  
  # index k = h(h-1)/2, with h = floor(n/2) + 1
  h <- n %/% 2L + 1L
  k_idx <- h * (h - 1L) / 2L   # 1-based
  base <- d[k_idx]
  
  # asymptotic consistency constant for sigma
  scale <- 2.21914446598508
  
  # Qn0 = scale * base (no finite-sample correction)
  scale * base
}

Sn_raw <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (n < 2L) return(0)
  
  row_meds <- numeric(n)
  for (i in 1:n) {
    d <- abs(x[i] - x)
    d <- sort(d)
    if (n %% 2L == 1L) {
      # odd n
      m <- d[(n + 1L) / 2L]
    } else {
      # even n, stable median
      a <- d[n / 2L]
      b <- d[n / 2L + 1L]
      m <- a + (b - a) / 2
    }
    row_meds[i] <- m
  }
  
  row_meds <- sort(row_meds)
  if (n %% 2L == 1L) {
    base <- row_meds[(n + 1L) / 2L]
  } else {
    a <- row_meds[n / 2L]
    b <- row_meds[n / 2L + 1L]
    base <- a + (b - a) / 2
  }
  
  # asymptotic consistency constant
  scale <- 1.19259855312321
  
  # Sn0 = scale * base (no finite-sample correction)
  scale * base
}

## ---------- Monte Carlo estimation for a vector of n ----------

estimate_dn_cn <- function(
    n_values,
    B = 100000,
    estimate_Qn = TRUE,
    estimate_Sn = TRUE
) {
  n_values <- sort(unique(as.integer(n_values)))
  res_Q <- NULL
  res_S <- NULL
  
  if (estimate_Qn) {
    dn_vals <- numeric(length(n_values))
  }
  if (estimate_Sn) {
    cn_vals <- numeric(length(n_values))
  }
  
  for (idx in seq_along(n_values)) {
    n <- n_values[idx]
    cat("n =", n, "(B =", B, ")\n")
    
    if (estimate_Qn) {
      sumQ <- 0
    }
    if (estimate_Sn) {
      sumS <- 0
    }
    
    for (b in 1:B) {
      x <- rnorm(n)
      
      if (estimate_Qn) {
        sumQ <- sumQ + Qn_raw(x)
      }
      if (estimate_Sn) {
        sumS <- sumS + Sn_raw(x)
      }
    }
    
    if (estimate_Qn) {
      meanQ0 <- sumQ / B
      dn_vals[idx] <- 1 / meanQ0
    }
    if (estimate_Sn) {
      meanS0 <- sumS / B
      cn_vals[idx] <- 1 / meanS0
    }
  }
  
  if (estimate_Qn) {
    res_Q <- data.frame(
      n = n_values,
      dn_hat = dn_vals
    )
  }
  if (estimate_Sn) {
    res_S <- data.frame(
      n = n_values,
      cn_hat = cn_vals
    )
  }
  
  list(
    Qn = res_Q,
    Sn = res_S
  )
}

## ---------- Example usage and file output ----------

## Choose n values
n_values <- c(2:12, 15, 20, 30, 50)

## Number of Monte Carlo replications per n
B <-800000  # increase for higher precision if needed

res <- estimate_dn_cn(
  n_values = n_values,
  B = B,
  estimate_Qn = TRUE,
  estimate_Sn = TRUE
)

cat("\nEstimated dn for Qn:\n")
print(res$Qn)

cat("\nEstimated cn for Sn:\n")
print(res$Sn)

## ---------- Save factors to files in the script directory ----------

qn_path <- file.path(script_dir, "qn_dn_factors.txt")
sn_path <- file.path(script_dir, "sn_cn_factors.txt")

if (!is.null(res$Qn)) {
  lines_qn <- sprintf("%d %.10f", res$Qn$n, res$Qn$dn_hat)
  writeLines(lines_qn, con = qn_path)
}

if (!is.null(res$Sn)) {
  lines_sn <- sprintf("%d %.10f", res$Sn$n, res$Sn$cn_hat)
  writeLines(lines_sn, con = sn_path)
}

cat("\nWritten:\n  ", qn_path, "\n  ", sn_path, "\n")
