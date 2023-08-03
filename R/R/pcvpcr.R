#' Procrustes cross-validation for PCR models
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
#' @param cv.scope
#' scope for center/scale operations inside CV loop: 'global' — using globally computed mean and std
#' or 'local' — recompute new for each local calibration set.
#'
#' @details
#' The method computes pseudo-validation matrix Xpv, based on PCR decomposition of calibration
#' set {X, y} and cross-validation.
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
#' Parameter `cv.scope` influences how the Procrustean rule is met. In case of "global" scope,
#' the rule will be met strictly - error of predictions for PV-set and the global model will be
#' identical to the error from conventional cross-validation. In case of "local" scope, every
#' local model will have its own center and hence the rule will be almost met (the errors will
#' be close but not identical).
#'
#' @return
#' Pseudo-validation matrix (same size as X) with an additional attribute, `D` which contains the
#' scaling coefficients (ck/c)
#'
#' @references
#' 1. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev. Procrustes cross-validation of multivariate
#' regression models. Analytica Chimica Acta, 1255 (2022)
#' [https://doi.org/10.1016/j.aca.2023.341096]
#'
#' @examples
#'
#' # load NIR spectra of Corn samples
#' data(corn)
#' X <- corn$spectra
#' Y <- corn$moisture
#'
#' # generate Xpv set based on PCA decomposition with A = 20 and venetian blinds split with 4 segments
#' Xpv <- pcvpcr(X, Y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))
#'
#' # show the original spectra and the PV-set (as is and mean centered)
#' oldpar <- par(mfrow = c(2, 2))
#' matplot(t(X), type = "l", lty = 1, main = "Original data")
#' matplot(t(Xpv), type = "l", lty = 1, main = "PV-set")
#' matplot(t(scale(X, scale = FALSE)), type = "l", lty = 1, main = "Original data (mean centered)")
#' matplot(t(scale(Xpv, scale = FALSE)), type = "l", lty = 1, main = "PV-set (mean centered)")
#' par(oldpar)
#'
#' @export
pcvpcr <- function(X, Y, ncomp = min(nrow(X) - 1, ncol(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE, cv.scope = "global") {

   funlist <- list(

      # computes global PCR model
      getglobalmodel = function(X, Y, ncomp) {
         m <- svd(X, nv = ncomp, nu = ncomp)
         P <- m$v[, seq_len(ncomp), drop = FALSE]
         T <- X %*% P
         C <- t(solve(crossprod(T)) %*% crossprod(T, Y))
         I <- diag(1, ncol(X))
         Pi <- I - tcrossprod(P)
         return(list(P = P, C = C, Pi = Pi, ncomp = ncomp))
      },

      # computes local PCR model
      getlocalmodel = function(X.c, Y.c, m) {
         m.k <- svd(X.c, nu = ncomp, nv = ncomp)
         P.k <- m.k$v[, seq_len(ncomp), drop = FALSE]
         aa <- acos(colSums(m$P * P.k)) < (pi / 2)
         P.k <- P.k %*% diag(aa * 2 - 1, ncol(P.k), ncol(P.k))
         T.c <- X.c %*% P.k
         C.k <- t(solve(crossprod(T.c)) %*% crossprod(T.c, Y.c))
         return(list(P = P.k, C = C.k))
      },

      # computes PV-set for current segment
      getxpv = function(m, m.k, X.k) {

         # get global model parameters for current component and compute regression vector
         P <- m$P
         C <- as.numeric(m$C)

         # get local model parameters for current component and compute regression vector
         P.k <- m.k$P
         C.k <- as.numeric(m.k$C)

         # compute scores for PV subset
         T.k <- X.k %*% P.k
         D.k <- diag(C.k / C, length(C), length(C))

         T.pv <- T.k %*% D.k
         X.pv <- tcrossprod(T.pv, P)

         return(list(X = X.pv, D = D.k))
      },

      # computes vector with orthogonal distances
      getqk = function(X.k, m.k) {
         return( rowSums( (X.k - tcrossprod(X.k %*% m.k$P, m.k$P))^2 ))
      }
   )

   return(pcvreg(X, Y, ncomp, cv = cv, center = center, scale = scale, funlist = funlist, cv.scope = cv.scope))
}
