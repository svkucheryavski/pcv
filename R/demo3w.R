#####################################################################################
#                                                                                   #
#  This file contains demo code for Procrustes cross-validation of three-way data   #
#  in R. It uses simulated EEM (excitation-emission) data, generates PV-sets with   #
#  the three methods from the package and shows how to use a PV-set for validation  #
#  of a three-way model.                                                            #
#                                                                                   #
#  Check this for more details: https://github.com/svkucheryavski/pcv/              #
#                                                                                   #
#####################################################################################

rm(list = ls())
library(pcv)


#####################################################################################
# 0. Simulated dataset                                                              #
#####################################################################################

## The data are similar to fluorescence EEM spectra: each sample is a matrix with
## excitation wavelengths in rows and emission wavelengths in columns. Every sample
## is a mixture of three fluorophores with random concentrations. Each fluorophore has
## a Gaussian shaped excitation and emission profile. The peak positions shift slightly
## from sample to sample (so the data are not exactly trilinear, as real EEM data) and
## small random noise is added.

set.seed(42)
I <- 60                                    # number of samples
ex <- seq(250, 400, length.out = 31)       # excitation wavelengths (mode 2)
em <- seq(300, 500, length.out = 41)       # emission wavelengths (mode 3)

peak <- function(x, mu, sd) exp(-0.5 * ((x - mu) / sd)^2)
ex.pos <- c(280, 320, 360)                 # excitation peak positions of the fluorophores
em.pos <- c(340, 400, 450)                 # emission peak positions
ex.sd <- c(15, 20, 15)                     # excitation peak widths
em.sd <- c(20, 25, 20)                     # emission peak widths
A <- matrix(runif(I * 3, 0.2, 1), I, 3)    # concentrations

## three-way array with dimension I x J x K (samples x excitation x emission)
X <- array(0, c(I, length(ex), length(em)))
for (i in seq_len(I)) {
   for (f in 1:3) {
      b <- peak(ex, ex.pos[f] + rnorm(1, sd = 4), ex.sd[f])
      c <- peak(em, em.pos[f] + rnorm(1, sd = 4), em.sd[f])
      X[i, , ] <- X[i, , ] + A[i, f] * outer(b, c)
   }
}
X <- X + array(rnorm(length(X), sd = 0.01), dim(X))


#####################################################################################
# 1. Generate PV-sets                                                               #
#####################################################################################

## All three methods take the calibration set as three-way array and return PV-set
## of the same size, in the original units. Parameter cv defines the split into segments
## in the same way as in pcvpca(). Parameter seed makes the result reproducible.

## FU (full unfolding): every sample is unfolded to a vector of J x K values and PV-set
## is generated with PCA based PCV. The number of components must be large enough,
## so the PCA subspace contains the subspace of the three-way models you validate.
Xpv.fu <- pcv3wfu(X, ncomp = 20, cv = list("ven", 4), seed = 1)

## BU (block unfolding): the array is split into J matrices I x K (mode = 2) or
## K matrices I x J (mode = 3), and PV-set is generated for each matrix separately
## with the same split of the samples.
Xpv.bu <- pcv3wbu(X, mode = 2, ncomp = 10, cv = list("ven", 4), seed = 1)

## HSV (HO-SVD based): the array is modeled with truncated HO-SVD in modes 2 and 3,
## with ncomp = c(R2, R3) components. Use more components than the largest model you
## plan to validate. If the last components are close in size, their order can change
## when some samples are removed, and then the PV-set reproduces cross-validation
## only for some smaller numbers of components. In this case you get a warning.
Xpv.hsv <- pcv3whsv(X, ncomp = c(6, 6), cv = list("ven", 4), seed = 1)

## HSV PV-set has attributes, e.g. the number of components actually used and the
## numbers of components for which the PV-set reproduces cross-validation exactly
attr(Xpv.hsv, "ncomp")
attr(Xpv.hsv, "valid.subranks")

## show the first sample of calibration set and of the PV-sets as heatmaps
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
showeem <- function(S, main) {
   image(ex, em, S, col = hcl.colors(64, "viridis"), xlab = "Excitation, nm",
      ylab = "Emission, nm", main = main)
}
showeem(X[1, , ], "Calibration set")
showeem(Xpv.fu[1, , ], "PV-set (FU)")
showeem(Xpv.bu[1, , ], "PV-set (BU, mode 2)")
showeem(Xpv.hsv[1, , ], "PV-set (HSV)")


#####################################################################################
# 2. Use HSV PV-set for validation of HO-SVD model                                  #
#####################################################################################

## PV-set is used in the same way as an independent validation set: you make a model
## for the calibration set, project the PV-set onto this model and compare the results.
## Here we use the HO-SVD (Tucker2) model, the model behind the HSV method, written
## in a few lines below. For PARAFAC or Tucker3 models you can use e.g. package
## "multiway" in the same way: fit the model on X, apply it to Xpv.

## fit HO-SVD model with R2 x R3 components to array X (with or without mean centering)
hosvd <- function(X, R2, R3, center = TRUE) {
   M <- if (center) apply(X, c(2, 3), mean) else matrix(0, dim(X)[2], dim(X)[3])
   Xc <- sweep(X, c(2, 3), M, "-")
   d <- dim(X)
   U2 <- svd(matrix(aperm(Xc, c(2, 3, 1)), d[2]), nu = R2, nv = 0)$u
   U3 <- svd(matrix(aperm(Xc, c(3, 2, 1)), d[3]), nu = R3, nv = 0)$u
   m <- list(M = M, U2 = U2, U3 = U3)

   # SVD of the calibration scores is needed to compute score distances
   T <- t(sapply(seq_len(d[1]), function(i) as.vector(crossprod(U2, Xc[i, , ]) %*% U3)))
   s <- svd(T, nu = 0)
   m$V <- s$v
   m$s <- s$d
   return(m)
}

## project array X onto the model, compute score distance (h) and orthogonal distance (q)
hosvd.predict <- function(m, X) {
   Xc <- sweep(X, c(2, 3), m$M, "-")
   res <- t(sapply(seq_len(dim(X)[1]), function(i) {
      Xi <- Xc[i, , ]
      Ti <- crossprod(m$U2, Xi) %*% m$U3
      h <- sum((crossprod(m$V, as.vector(Ti)) / m$s)^2)
      q <- sum(Xi^2) - sum(Ti^2)
      c(h = h, q = q)
   }))
   as.data.frame(res)
}

## distance plot: calibration set and PV-set, both distances scaled by the mean values
## for the calibration set, similar to the distance plots in DD-SIMCA
plotdist <- function(m, X, Xpv, main) {
   rcal <- hosvd.predict(m, X)
   rpv <- hosvd.predict(m, Xpv)
   h0 <- mean(rcal$h)
   q0 <- mean(rcal$q)
   xlim <- c(0, max(rcal$h / h0, rpv$h / h0))
   ylim <- c(0, max(rcal$q / q0, rpv$q / q0))
   plot(rcal$h / h0, rcal$q / q0, pch = 22, col = "gray", bg = "gray", xlim = xlim, ylim = ylim,
      xlab = "Score distance, h/h0", ylab = "Orthogonal distance, q/q0", main = main)
   points(rpv$h / h0, rpv$q / q0, pch = 4, col = "red")
   legend("topright", c("calibration", "PV-set"), pch = c(22, 4), col = c("gray", "red"),
      pt.bg = c("gray", NA), bty = "n")
}

## Two models are shown: with 2 x 2 components and with 6 x 6 components. For the
## first model the distances of the PV-set are similar to the calibration set. For
## the second model the score distances of the PV-set are much larger, as they are
## in cross-validation: the extra components describe variation which changes when
## some of the samples are removed, so the model is too complex.
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plotdist(hosvd(X, 2, 2), X, Xpv.hsv, "HO-SVD model, 2 x 2 components")
plotdist(hosvd(X, 6, 6), X, Xpv.hsv, "HO-SVD model, 6 x 6 components")


#####################################################################################
# 3. PV-set reproduces cross-validation                                             #
#####################################################################################

## Projected onto the global model, HSV PV-set gives the same scores and orthogonal
## distances as cross-validation of the HO-SVD model. Score distances are also the same
## if cross-validation computes them with the score metric of the global model (after
## matching local components to the global ones), but not with the metric of each local
## model. Below we run the cross-validation directly and compare the orthogonal
## distances for the 2 x 2 model.

cvind <- attr(Xpv.hsv, "cvind")
Xc <- sweep(X, c(2, 3), apply(X, c(2, 3), mean), "-")   # global centering, as in pcv3whsv()
q.cv <- numeric(I)
for (k in unique(cvind)) {
   mk <- hosvd(Xc[cvind != k, , , drop = FALSE], 2, 2, center = FALSE)
   q.cv[cvind == k] <- hosvd.predict(mk, Xc[cvind == k, , , drop = FALSE])$q
}
q.pv <- hosvd.predict(hosvd(X, 2, 2), Xpv.hsv)$q

par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(q.cv, q.pv, xlab = "Orthogonal distance, cross-validation",
   ylab = "Orthogonal distance, PV-set", main = "HSV PV-set vs. cross-validation (2 x 2)")
abline(0, 1, col = "gray", lty = 2)
max(abs(q.pv - q.cv))
