#!/usr/bin/env Rscript
# Bedarf: install.packages("lmtest")
library(lmtest)

# Daten einlesen
y <- scan("Documents/Develop_current/personal/SwiftyStatsGH/R-Scripts/lew.txt", sep = ",")

# Modell wählen: Intercept + Trend (passt zum Swift-Beispiel)
t <- seq_along(y)
fit <- lm(y ~ t)

# Durbin–Watson-Test (dwtest)
dw_greater <- dwtest(fit, alternative = "greater")
dw_less    <- dwtest(fit, alternative = "less")
dw_two     <- dwtest(fit, alternative = "two.sided")

cat(sprintf("DW statistic: %.4f\n", dw_two$statistic))
cat(sprintf("p(rho>0): %.6f\n", dw_greater$p.value))
cat(sprintf("p(rho<0): %.6f\n", dw_less$p.value))
cat(sprintf("p(two-sided): %.6f\n", dw_two$p.value))