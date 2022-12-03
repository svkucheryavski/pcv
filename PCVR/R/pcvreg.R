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
#' list with functions for particular implementatio
#'
#' @description
#' This is a generic method, use `pcvpls()` or `pcvpcr()` instead.
#'
#' @importFrom stats sd
pcvreg <- function(X, Y, ncomp = min(nrow(X) - 1, ncol(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE, funlist = list()) {

   if (is.null(dim(Y))) {
      dim(Y) <- c(nrow(X), 1)
   }

   # keep names if any
   attrs <- attributes(X)

   # compute center and scale values for predictors
   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale) apply(X, 2, sd) else rep(1, ncol(X))

   # compute center and scale values for responses
   mY <- if (center) apply(Y, 2, mean) else rep(0, ncol(Y))
   sY <- if (scale) apply(Y, 2, sd) else rep(1, ncol(Y))

   # autoscale both
   X <- scale(X, center = mX, scale = sX)
   Y <- scale(Y, center = mY, scale = sY)

   # get indices for cross-validation
   ind <- pcvcrossval(cv, nrow(X), Y[, 1])

   # get number of segments
   nseg <- max(ind)

   # correct maximum number of components
   max.ncomp <- min(nrow(X) - ceiling(nrow(X) / nseg) - 1, ncol(X))
   if (ncomp > max.ncomp) {
      ncomp <- max.ncomp
   }

   # prepare empty matrix for pseudo-validation set
   Xpv <- matrix(0, nrow(X), ncol(X))

   # compute global model
   m <- funlist$getglobalmodel(X, Y, ncomp)
   D <- matrix(0, nseg, ncomp)

   # loop for computing the PV set
   for (k in seq_len(nseg)) {

      # split data to calibration and validation
      ind.c <- ind != k
      ind.k <- ind == k

      X.c <- X[ ind.c, , drop = FALSE]
      X.k <- X[ ind.k, , drop = FALSE]
      Y.c <- Y[ ind.c, , drop = FALSE]

      # compute local model for current segment
      m.k <- funlist$getlocalmodel(X.c, Y.c, m)

      # compute explained part of PV-set for current segment
      pvres <- funlist$getxpv(m, m.k, X.k)
      Xpv.hat <- pvres$X
      D[k, ] <- diag(pvres$D)

      # compute the orthogonal part of PV-set
      Xpv.orth <- getxpvorth(funlist$getqk(X.k, m.k), X.k, m$Pi)

      # add the orthogonal part
      Xpv[ind.k, ] <- Xpv.hat + Xpv.orth
   }

   # uncenter and unscale the data
   Xpv <- sweep(Xpv, 2, sX, "*")
   Xpv <- sweep(Xpv, 2, mX, "+")

   # add initial attributes
   attributes(Xpv) <- attrs

   # add additional attribute with diagonal elements
   attr(Xpv, "D") <- D

   return(Xpv)
}

