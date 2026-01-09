#!/usr/bin/env Rscript

# Scheffe post-hoc test for the provided groups (one-way ANOVA), using DescTools.
library("DescTools")
g1 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g2 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g3 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g4 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g5 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g6 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g7 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g8 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g9 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
g10 <- c(1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03)
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
