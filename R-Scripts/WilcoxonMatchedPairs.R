## Example data from ?wilcox.test
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

## One-sided exact Wilcoxon signed-rank test
(wt <- wilcoxsign_test(y1 ~ y2, distribution = "exact",
                       alternative = "two.sided"))
statistic(wt, type = "linear")
midpvalue(wt)
pvalue(wt)

## Example data from ?wilcox.test
y1 <- c(63,65,71,75,90,75,68,74,62,73)
y2 <- c(80,78,96,87,88,96,82,73,77,79)

## One-sided exact Wilcoxon signed-rank test
(wt <- wilcoxsign_test(y1 ~ y2, distribution = "exact",
                       alternative = "two.sided", correction = TRUE, zero.method = "Pratt"))
statistic(wt, type = "linear")
midpvalue(wt)
pvalue(wt)
