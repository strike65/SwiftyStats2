#!/usr/bin/env Rscript

# Scheffe post-hoc test for the provided groups (one-way ANOVA), using DescTools.
library("DescTools")
g1 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g2 <- c(1.20, 1.18, 1.25, 1.14, 1.22, 1.19, 1.28, 1.16, 1.24, 1.21)
g3 <- c(1.35, 1.31, 1.29, 1.38, 1.33, 1.36, 1.40, 1.27, 1.34, 1.32)
g4 <- c(1.45, 1.50, 1.47, 1.42, 1.55, 1.49, 1.53, 1.44, 1.51, 1.46)
g5 <- c(1.60, 1.58, 1.63, 1.66, 1.55, 1.61, 1.64, 1.57, 1.62, 1.59)
g6 <- c(1.72, 1.68, 1.75, 1.70, 1.73, 1.67, 1.78, 1.71, 1.74, 1.69)
g7 <- c(1.85, 1.80, 1.88, 1.83, 1.90, 1.82, 1.87, 1.81, 1.89, 1.84)
g8 <- c(2.00, 1.95, 2.05, 1.98, 2.02, 1.97, 2.08, 1.96, 2.04, 1.99)
g9 <- c(2.15, 2.10, 2.18, 2.12, 2.20, 2.11, 2.17, 2.09, 2.19, 2.14)
g10 <- c(2.30, 2.25, 2.35, 2.28, 2.32, 2.27, 2.38, 2.26, 2.34, 2.29)
g11 <- c(2.45, 2.40, 2.50, 2.43, 2.48, 2.41, 2.55, 2.42, 2.49, 2.44)

values <- c(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11)
groups <- factor(rep(
  c("g01","g02","g03","g04","g05","g06","g07","g08","g09","g10","g11"),
  times = c(length(g1), length(g2), length(g3), length(g4), length(g5),
            length(g6), length(g7), length(g8), length(g9), length(g10), length(g11))
))
df <- data.frame(value = values, group = groups)

anova_fit <- aov(value ~ group, data = df)
print(summary(anova_fit))

# Ensure DescTools is available
if (!requireNamespace("DescTools", quietly = TRUE)) {
  install.packages("DescTools", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(DescTools))

# Scheffe post-hoc
# By default: all pairwise comparisons; set conf.level as needed.
scheffe_result <- ScheffeTest(anova_fit, which = "group", conf.level = 0.95)
print(scheffe_result)

# Optional: unified PostHocTest interface (same result family)
# posthoc_result <- PostHocTest(anova_fit, which = "group", method = "scheffe", conf.level = 0.95)
# print(posthoc_result)
