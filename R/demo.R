# This is demo code for the new version of the Procrustes cross-validation method (2022)
# It requires "pcv" and "mdatools" packages, if you do not have one
# Run: install.packages("pcv", "mdatools")
#
# The code creates plots similar to what you can see on Figures with Corn example
# in the paper describing the new version, submitted in 2022. The Corn dataset was
# downloaded from here: https://eigenvector.com/resources/data-sets/
#

rm(list = ls())
library(mdatools)
library(pcv)

# load original spectra and response values for Corn dataset
data(corn)
X = corn$spectra
Y = corn$moisture

# 1. PCA based examples

## create pseudo-validation set
Xpv <- pcvpca(X, ncomp = 30, cv = list("ven", 4))

## show plot with original and generated spectra
par(mfrow = c(2, 2))
mdaplot(X, type = "l", cgroup = as.numeric(Y), main = "Original")
mdaplot(Xpv, type = "l", cgroup = as.numeric(Y), main = "Pseudo-validation")
mdaplot(prep.autoscale(X), type = "l", cgroup = as.numeric(Y), main = "Original (mean centered)")
mdaplot(prep.autoscale(Xpv), type = "l", cgroup = as.numeric(Y), main = "Pseudo-validation (mean centered)")

## make PCA model for calibration set
m <- pca(X, 20)

## project pseudo-validation set to the model
res_pv <- predict(m, Xpv)

## show Distance plot for A = 2 and A = 20
par(mfrow = c(1, 2))
plotResiduals(m, res = list(cal = m$res$cal, pv = res_pv), pch = c(22, 4), ylim = c(0, 12), ncomp = 2, main = "Distance plot (A = 2)")
plotResiduals(m, res = list(cal = m$res$cal, pv = res_pv), pch = c(22, 4), ylim = c(0, 12), ncomp = 20, main = "Distance plot (A = 20)")


# 2. PLS based examples

## create pseudo-validation set for PLS model
Xpv <- pcvpls(X, Y, ncomp = 30, cv = list("ven", 4))

## create PLS model using xpv as validation set
m <- pls(X, Y, ncomp = 20, x.test = Xpv, y.test = Y)

## show the main plots
par(mfrow = c(2, 2))
plotXCumVariance(m)
plotYCumVariance(m)
plotRMSE(m)
plotPredictions(m, 10)

## show heatmap for elements of D - it is returned as attribute of PC-set
## this plot must be identical to the first plot in the figure with heatmaps shown in the manuscript

D <- attr(Xpv, "D")

colmap <- mdaplot.getColors(256, colmap = c("blue", "white", "red"))
par(mfrow = c(2, 1))
plotD(Xpv)

boxplot(D, range = 100, col = "gray", border = "gray", ylim = c(0, 2), xlab = "Components", ylab = expression(c[k]/c))
abline(h = 1, col = "black", lty = 2)
par(mar = c(5, 4, 2, 2))
