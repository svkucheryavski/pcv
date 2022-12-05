#' Procrustes cross-validation for PCA models
#'
#' @param X
#' matrix with predictors from the training set.
#' @param ncomp
#' number of components to use (more than the expected optimal number).
#' @param center
#' logical, to center or not the data sets
#' @param scale
#' logical, to scale or not the data sets
#' @param cv
#' which split method to use for cross-validation (see description for details).
#'
#' @details
#' The method computes pseudo-validation matrix Xpv, based on PCA decomposition of calibration
#' set `X` and cross-validation. See description of the method in [1].
#'
#' Parameter `cv` defines how to split the rows of the training set. The split is similar
#' to cross-validation splits, as PCV is based on cross-validation. This parameter can have
#' the following values:
#'
#' 1. A number (e.g. `cv = 4`). In this case this number specifies number of segments for random
#' splits, except `cv = 1` which is a special case for leave-one-out (full cross-validation).
#'
#' 2. A list with 2 values: `list("name", nseg)`. In this case `"name"` defines the way to make
#' the split, you can select one of the following: `"loo"` for leave-one-out, `"rand"` for random
#' splits or `"ven"` for Venetian blinds (systematic) splits. The second parameter, `nseg`, is a
#' number of segments for splitting the rows into. For example, `cv = list("ven", 4)`, shown in
#' the code examples above, tells PCV to use Venetian blinds splits with 4 segments.
#'
#' 3. A vector with integer numbers, e.g. `cv = c(1, 2, 3, 1, 2, 3, 1, 2, 3)`. In this case number
#' of values in this vector must be the same as number of rows in the training set. The values
#' specify which segment a particular row will belong to. In case of the example shown here, it
#' is assumed that you have 9 rows in the calibration set, which will be split into 3 segments.
#' The first segment will consist of measurements from rows 1, 4 and 7.
#'
#' @return
#' Matrix with PV-set (same size as X)
#'
#' @references
#' 1. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev. Procrustes cross-validation of multivariate
#' regression models. Submitted, 2022.
#'
#' @examples
#'
#' # load NIR spectra of Corn samples
#' data(corn)
#' X <- corn$spectra
#'
#' # generate Xpv set based on PCA decomposition with A = 20 and venetian blinds split with 4 segments
#' Xpv <- pcvpca(X, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))
#'
#' # show the original spectra and the PV-set (as is and mean centered)
#' oldpar <- par(mfrow = c(2, 2))
#' matplot(t(X), type = "l", lty = 1, main = "Original data")
#' matplot(t(Xpv), type = "l", lty = 1, main = "PV-set")
#' matplot(t(scale(X, scale = FALSE)), type = "l", lty = 1, main = "Original data (mean centered)")
#' matplot(t(scale(Xpv, scale = FALSE)), type = "l", lty = 1, main = "PV-set (mean centered)")
#' par(oldpar)
#'
#' @importFrom stats sd
#'
#' @export
pcvpca <- function(X, ncomp = min(nrow(X) - 1, col(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE) {

   # keep names if any
   attrs <- attributes(X)

   mX <- apply(X, 2, mean)
   sX <- if (scale) apply(X, 2, sd) else rep(1, ncol(X))

   # autoscale the calibration set
   X <- scale(X, center = mX, scale = sX)

   # indices for cross-validation
   ind <- pcvcrossval(cv, nrow(X), seq_len(nrow(X)))

   # get number of segments
   nseg <- max(ind)

   # correct maximum number of components
   max.ncomp <- min(nrow(X) - ceiling(nrow(X) / nseg) - 1 , ncol(X))
   if (ncomp > max.ncomp) {
      ncomp <- max.ncomp
   }

   # create a global model
   P <- svd(X)$v[, seq_len(ncomp), drop = FALSE]
   Pi <- diag(1, nrow(P)) - tcrossprod(P)

   # prepare empty matrix for pseudo-validation set
   Xpv <- matrix(0, nrow(X), ncol(X))

   # loop for computing the PV set
   for (k in seq_len(nseg)) {

      # split data to calibration and validation
      ind.c <- ind != k
      ind.k <- ind == k

      X.c <- X[ ind.c, , drop = FALSE]
      X.k <- X[ ind.k, , drop = FALSE]

      # get loadings for local model and rotation matrix between global and local models
      P.k <- svd(X.c, nv = ncomp)$v[, seq_len(ncomp), drop = FALSE]

      # correct direction of loadings for local model
      a <- acos(colSums(P * P.k)) < pi / 2
      P.k <- P.k %*% diag(a * 2 - 1, ncol(P), ncol(P))

      # get scores and residuals by projection local validation set to the local model
      T.k <- X.k %*% P.k
      E.k <- X.k - tcrossprod(T.k, P.k)
      q.k <- rowSums(E.k^2)

      # compute the parallel component of Xpv
      Xpv.hat <- tcrossprod(T.k, P)

      # compute the orthogonal component of Xpv
      Xpv.orth <- getxpvorth(q.k, X.k, Pi)

      # create and save the Xpv
      Xpv[ind.k, ] <- Xpv.hat + Xpv.orth
   }

   # uscenter and unscale the data
   Xpv <- sweep(Xpv, 2, sX, "*")
   Xpv <- sweep(Xpv, 2, mX, "+")

   attributes(Xpv) <- attrs
   return(Xpv)
}


