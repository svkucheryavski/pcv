#' Procrustes cross-validation for multivariate regression models
#'
#' @param X
#' matrix with predictors from the training set.
#' @param Y
#' vector with response values from the training set.
#' @param ncomp
#' number of components to use (more than the expected optimal number).
#' @param center
#' logical, to center or not the data sets
#' @param scale
#' logical, to scale or not the data sets
#' @param cv
#' which split method to use for cross-validation (see description of method `pcvpls()` for details).
#' @param funlist
#' list with functions for particular implementation
#' @param cv.scope
#' scope for center/scale operations inside CV loop: 'global' — using globally computed mean and std
#' or 'local' — recompute new for each local calibration set.
#'
#' @description
#' This is a generic method, use `pcvpls()` or `pcvpcr()` instead.
#'
#' @importFrom stats sd
pcvreg <- function(X, Y, ncomp = min(nrow(X) - 1, ncol(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE, funlist = list(), cv.scope = "global") {

   nRows <- nrow(X)
   nPred <- ncol(X)

   if (is.null(dim(Y))) {
      dim(Y) <- c(nRows, 1)
   }

   nResp <- ncol(Y)

   # keep names if any
   attrs <- attributes(X)

  # get indices for cross-validation
   ind <- pcvcrossval(cv, nRows, Y[, 1])

   # get number of segments
   nSeg <- max(ind)

   # correct maximum number of components
   max.ncomp <- min(nrow(X) - ceiling(nRows / nSeg) - 1, nPred)
   if (ncomp > max.ncomp) {
      ncomp <- max.ncomp
   }

   # compute global center and scale values for predictors
   mX <- rep(0, nPred)
   mY <- rep(0, nResp)
   if (center) {
      mX <- apply(X, 2, mean)
      mY <- apply(Y, 2, mean)
   }

   sX <- rep(1, nPred)
   sY <- rep(1, nResp)
   if (scale) {
      sX <- apply(X, 2, sd)
      sY <- apply(Y, 2, sd)
   }

   # autoscale globally
   X <- scale(X, center = mX, scale = sX)
   Y <- scale(Y, center = mY, scale = sY)
   m <- funlist$getglobalmodel(X, Y, ncomp)

   # prepare empty matrix for pseudo-validation set
   Xpv <- matrix(0, nRows, nPred)
   D <- matrix(0, nSeg, ncomp)

   # this will be used only for PLS case
   R <- array(0, dim = c(nSeg, nPred, ncomp))

   # loop for computing the PV set
   for (k in seq_len(nSeg)) {

      # split data to calibration and validation
      ind.c <- ind != k
      ind.k <- ind == k

      X.k <- X[ind.k, , drop = FALSE]
      X.c <- X[ind.c, , drop = FALSE]
      Y.c <- Y[ind.c, , drop = FALSE]

      # if cv.scope is local autoscale locally
      if (cv.scope == "local") {

         mXl <- rep(0, nPred)
         mYl <- rep(0, nResp)
         if (center) {
            mXl <- apply(X.c, 2, mean)
            mYl <- apply(Y.c, 2, mean)
         }

         sXl <- rep(1, nPred)
         sYl <- rep(1, nResp)
         if (scale) {
            sXl <- apply(X.c, 2, sd)
            sYl <- apply(Y.c, 2, sd)
         }

         X.c <- scale(X.c, center = mXl, scale = sXl)
         Y.c <- scale(Y.c, center = mYl, scale = sYl)
         X.k <- scale(X.k, center = mXl, scale = sXl)
      }

      # compute local model for current segment
      m.k <- funlist$getlocalmodel(X.c, Y.c, m)

      # compute explained part of PV-set for current segment
      pvres <- funlist$getxpv(m, m.k, X.k)
      Xpv.hat <- pvres$X
      D[k, ] <- diag(pvres$D)

      # in case of PLS regression we save R as well
      if (!is.null(pvres$R)) R[k, , ] <- pvres$R

      # compute the orthogonal part of PV-set
      Xpv.orth <- getxpvorth(funlist$getqk(X.k, m.k), X.k, m$Pi)

      # add the orthogonal part
      Xpv[ind.k, ] <- Xpv.hat + Xpv.orth
   }

   # uncenter and unscale the data for global scope
   Xpv <- sweep(Xpv, 2, sX, "*")
   Xpv <- sweep(Xpv, 2, mX, "+")

   # add initial attributes
   attributes(Xpv) <- attrs

   # add additional attribute with diagonal elements
   attr(Xpv, "D") <- D

   # in case of PLS regression we save R as well
   if (!is.null(pvres$R)) {
      attr(Xpv, "Rk") <- R
      attr(Xpv, "R") <- m$R
   }

   return(Xpv)
}

