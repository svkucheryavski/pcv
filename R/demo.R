# This is demo code for Procrustes cross-validation method
# It requires "mdatools" package, if you do not have one
# Run: install.packages("mdatools")

library(mdatools)
rm(list = ls())
source("pcv.R")

# load original spectra
load("nirsim.RData")
write.table(spectra, file = "nirsim.csv", sep = ";", dec = ",")

# create pseudo-validation set
spectra_pv <- pcv(spectra, ncomp = 6, nseg = 4)

# show plot with original and generated spectra
par(mfrow = c(2, 1))
mdaplot(spectra, type = "l", main = "Original")
mdaplot(spectra_pv, type = "l", main = "Pseudo-validation")

# make PCA model for calibration set
m <- pca(spectra, 6)

# project pseudo-validation set to the model
res_pv <- predict(m, spectra_pv)

# show extreme plot for A = 4
par(mfrow = c(1, 2))
plotExtreme(m, res = m$res$cal, comp = 4, main = "Extreme (cal)")
plotExtreme(m, res = res_pv, comp = 4, main = "Extreme (pv)")

