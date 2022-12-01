# This is demo code for the new version of the Procrustes cross-validation method (2022)
# It requires "mdatools" package, if you do not have one
# Run: install.packages("mdatools")
#
# The code creates plots similar to what you can see on Figures with Corn example
# in the paper describing the new version, submitted in 2022. The Corn dataset was
# downloaded from here: https://eigenvector.com/resources/data-sets/
#

rm(list = ls())
library(mdatools)
source("pcv.R")

# load original spectra and response values for Corn dataset
load("corn.RData")

# do SNV preprocessing
x = prep.snv(x)

# 1. PCA based examples

## create pseudo-validation set
xpv <- pcvpca(x, ncomp = 30, cv = list("ven", 4))

## show plot with original and generated spectra
par(mfrow = c(2, 2))
mdaplot(x, type = "l", cgroup = as.numeric(y), main = "Original")
mdaplot(xpv, type = "l", cgroup = as.numeric(y), main = "Pseudo-validation")
mdaplot(prep.autoscale(x), type = "l", cgroup = as.numeric(y), main = "Original (mean centered)")
mdaplot(prep.autoscale(xpv), type = "l", cgroup = as.numeric(y), main = "Pseudo-validation (mean centered)")

## make PCA model for calibration set
m <- pca(x, 20)

## project pseudo-validation set to the model
res_pv <- predict(m, xpv)

## show Distance plot for A = 2 and A = 20
par(mfrow = c(1, 2))
plotResiduals(m, res = list(cal = m$res$cal, pv = res_pv), pch = c(22, 4), ylim = c(0, 12), ncomp = 2, main = "Distance plot (A = 2)")
plotResiduals(m, res = list(cal = m$res$cal, pv = res_pv), pch = c(22, 4), ylim = c(0, 12), ncomp = 20, main = "Distance plot (A = 20)")


# 2. PLS based examples

## create pseudo-validation set for PLS model
xpv <- pcvpls(x, y, ncomp = 30, cv = list("ven", 4))

## create PLS model using xpv as validation set
m <- pls(x, y, ncomp = 20, x.test = xpv, y.test = y)

## show the main plots
par(mfrow = c(2, 2))
plotXCumVariance(m)
plotYCumVariance(m)
plotRMSE(m)
plotPredictions(m, 10)

## show heatmap for elements of D - it is returned as attribute of PC-set
## this plot must be identical to the first plot in the figure with heatmaps shown in the manuscript

D <- attr(xpv, "D")

colmap <- mdaplot.getColors(256, colmap = c("blue", "white", "red"))
par(mfrow = c(2, 1))
par(mar = c(1, 4, 2, 1))
image(seq_len(30), seq_len(4), t(D), zlim = c(-2, 4), breaks = seq(-2, 4, length.out = 257),
   col = colmap, xlab = "Components (a)", ylab = "Segments (k)", xaxt = "n", main = "K = 4")
par(mar = c(5, 4, 1, 1))
boxplot(D, range = 100, col = "gray", border = "gray", ylim = c(0, 2), xlab = "Components", ylab = expression(c[k]/c))
abline(h = 1, col = "black", lty = 2)
par(mar = c(5, 4, 2, 2))
