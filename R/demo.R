# This is demo code for Procrustes cross-validation method
# It requires "mdatools" package, if you do not have one 
# Run: install.packages("mdatools") 

library(mdatools)
rm(list = ls())
source("pcv.R")

# load original spectra
load("nirsim.RData")

# create pseudo-validation set
spectra_pv <- pcv(spectra, ncomp = 8, nseg = 4)

# show plot with original and generated spectra
par(mfrow = c(2, 1))
matplot(spectra, type = "l", col = "blue", lty = 1, main = "Original")
matplot(spectra_pv, type = "l", col = "blue", lty = 1, main = "Pseudo-validation")

# make PCA model for calibration set
m <- pca(spectra, 6)

# project pseudo-validation set to the model
res_pv <- predict(m, spectra_pv)

# show distance plot for A = 2 and A = 5
par(mfrow = c(1, 2))
plotResiduals(m, ncomp = 2, res = list(cal = m$res$cal, pv = res_pv)
plotResiduals(m, ncomp = 5, res = list(cal = m$res$cal, pv = res_pv)

