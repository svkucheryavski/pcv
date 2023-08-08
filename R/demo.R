#####################################################################################
#                                                                                   #
#  This file contains demo code for implementation of Procrustes cross-validation   #
#  method in R. The code creates plots similar to what you can see on figures       #
#  with Corn example in this paper: https://doi.org/10.1016/j.aca.2023.341096       #
#                                                                                   #
#  Check this for more details: https://github.com/svkucheryavski/pcv/              #
#                                                                                   #
#####################################################################################

rm(list = ls())
library(mdatools)
library(pcv)

# load original spectra and response values for Corn dataset
data(corn)
X = corn$spectra
Y = corn$moisture

# 1. PCA based examples

## create pseudo-validation set using PCA based algorithm
Xpv <- pcvpca(X, ncomp = 30, cv = list("ven", 4))

## make PCA model for calibration set
m <- pca(X, 20)

## project calibration and PV-set to the PCA model
res_cal <- predict(m, X)
res_pv <- predict(m, Xpv)

## show distance plots for A = 2 and A = 20
## similar to ones shown in Supplementary materials
ylim <- c(0, 12)
pch <- c(22, 4)
par(mfrow = c(1, 2))
plotResiduals(m, res = list(cal = res_cal, pv = res_pv), pch = pch, ylim = ylim, ncomp = 2, main = "Distance plot (A = 2)")
plotResiduals(m, res = list(cal = res_cal, pv = res_pv), pch = pch, ylim = ylim, ncomp = 20, main = "Distance plot (A = 20)")


# 2. PCR based examples

## create pseudo-validation set using PCR based algorithm
Xpv <- pcvpcr(X, Y, ncomp = 20, cv = list("ven", 4))

## show plot with original and generated spectra
## similar to plots shown in Figure 2
par(mfrow = c(1, 2))
mdaplot(prep.autoscale(X), type = "l", cgroup = as.numeric(Y), main = "Calibration set (mean centered)")
mdaplot(prep.autoscale(Xpv), type = "l", cgroup = as.numeric(Y), main = "PV-set (mean centered)")


# 3. PLS based examples

## create pseudo-validation set for PLS model
Xpv <- pcvpls(X, Y, ncomp = 20, cv = list("ven", 4))

## create PLS model using xpv as validation set
m <- pls(X, Y, ncomp = 20, x.test = Xpv, y.test = Y)

## show the plots similar to Figure 5 in the manuscript
par(mfrow = c(1, 3))
plotRMSE(m)
plotPredictions(m, 10)
plotPredictions(m, 20)

## show heatmap for elements of D - it is returned as attribute of PC-set
## this plot must be identical to the first plot in the figure with heatmaps shown in the manuscript

D <- attr(Xpv, "D")

colmap <- mdaplot.getColors(256, colmap = c("blue", "white", "red"))
par(mfrow = c(2, 1))
plotD(Xpv)

boxplot(D, range = 100, col = "gray", border = "gray", ylim = c(0, 2), xlab = "Components", ylab = expression(c[k]/c))
abline(h = 1, col = "black", lty = 2)
