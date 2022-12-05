#' Procrustes cross-validation for PLS models
#'
#' @param X
#' matrix with predictors from the training set.
#' @param Y
#' vector or matrix with response values from the training set.
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
#' The method computes pseudo-validation matrix Xpv, based on PLS decomposition of calibration
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
#' @return
#' Pseudo-validation matrix (same size as X) with an additional attribute, `D` which contains the
#' scaling coefficients (ck/c)
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
#' Y <- corn$moisture
#'
#' # generate Xpv set based on PCA decomposition with A = 20 and venetian blinds split with 4 segments
#' Xpv <- pcvpls(X, Y, ncomp = 20, center = TRUE, scale = FALSE, cv = list("ven", 4))
#'
#' # show the original spectra and the PV-set (as is and mean centered)
#' oldpar <- par(mfrow = c(2, 2))
#' matplot(t(X), type = "l", lty = 1, main = "Original data")
#' matplot(t(Xpv), type = "l", lty = 1, main = "PV-set")
#' matplot(t(scale(X, scale = FALSE)), type = "l", lty = 1, main = "Original data (mean centered)")
#' matplot(t(scale(Xpv, scale = FALSE)), type = "l", lty = 1, main = "PV-set (mean centered)")
#' par(oldpar)
#'
#' # show the heatmap with the scaling coefficients
#' plotD(Xpv)
#'
#' @export
pcvpls <- function(X, Y, ncomp = min(nrow(X) - 1, ncol(X), 30), center = TRUE, scale = FALSE, cv = list("ven", 4)) {

   funlist <- list(

      # computes global PLS model
      getglobalmodel = function(X, Y, ncomp) {

         m <- simpls(X, Y, ncomp)
         P <- m$P
         C <- m$C
         R <- m$R
         Pi <- (diag(1, ncol(X)) - tcrossprod(R, P))
         return(list(P = P, Pi = Pi, C = C, R = R, ncomp = ncomp))
      },

      # computes local PLS model
      getlocalmodel = function(X.c, Y.c, m) {

         # local model
         m.k <- simpls(X.c, Y.c, m$ncomp)
         P.k <- m.k$P
         C.k <- m.k$C
         R.k <- m.k$R

         # correct angles for loadings and weights
         aa <- acos(colSums(normalize(m$R, 2) * normalize(R.k, 2))) < pi / 2
         R.k <- R.k %*% diag(aa * 2 - 1, ncol(R.k), ncol(R.k))
         P.k <- P.k %*% diag(aa * 2 - 1, ncol(P.k), ncol(P.k))
         C.k <- C.k %*% diag(aa * 2 - 1, ncol(C.k), ncol(C.k))

         return(list(P = P.k, R = R.k, C = C.k))
      },

      # computes PV-set for a segment k
      getxpv = function(m, m.k, X.k) {

         # get global model parameters for current component and compute regression vector
         R <- m$R
         C <- m$C
         P <- m$P

         # get local model parameters for current component and compute regression vector
         R.k <- m.k$R
         C.k <- m.k$C
         P.k <- m.k$P
         T.k <- X.k %*% R.k

         # compute the diagonal elements of the matrix D
         D.k <- diag(nrow = m$ncomp, ncol = m$ncomp)
         for (a in seq_len(m$ncomp)) {
            D.k[a, a] <- as.numeric(crossprod(C.k[, a], C[, a])) / as.numeric(crossprod(C[, a]))
         }

         T.pv <- T.k %*% D.k
         X.pv <- tcrossprod(T.pv, P)

         return(list(X = X.pv, D = D.k))
      },

      # computes the vector with orthogonal distances for local model
      getqk = function(X.k, m.k) {
         T.k <- X.k %*% m.k$R
         return( rowSums( (X.k - tcrossprod(T.k, m.k$P))^2 ))
      }
   )

   return(pcvreg(X, Y, ncomp, cv = cv, center = center, scale = scale, funlist = funlist))
}



#' Normalization rows or columns of a matrix
#' @param X
#' matrix with numeric values
#' @param dim
#' which dimension to normalize (1 for rows, 2 for columns)
#' @param weights
#' vector with normalization weights, by default 2-norm is used
#'
normalize <- function(X, dim = 1, weights = if(dim == 1) 1 / sqrt(rowSums(X^2)) else 1 / sqrt(colSums(X^2)) ) {
  sweep(X, dim, weights, FUN = "*")
}


#' SIMPLS algorithm
#'
#' @description
#' SIMPLS algorithm for calibration of PLS model
#'
#' @param X
#' a matrix with x values (predictors)
#' @param Y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#'
#' @return
#' a list with computed weights, x- and y-loadings for PLS regression model.
#'
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression.
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
#'
simpls <- function(X, Y, ncomp) {

   X <- as.matrix(X)
   Y <- as.matrix(Y)

   nobj  <- nrow(X)
   npred <- ncol(X)
   nresp <- ncol(Y)

   # initial estimation
   S <- crossprod(X, Y)
   M <- crossprod(X)

   # prepare space for results
   C <- matrix(0, nrow = nresp, ncol = ncomp)
   R <- V <- P <- matrix(0, nrow = npred, ncol = ncomp)
   TT <- U <- matrix(0, nrow = nobj, ncol = ncomp)


   # loop for each components
   for (a in seq_len(ncomp)) {

      r <- svd(S, nu = 1, nv = 0)$u
      t <- X %*% r

      tnorm <- sqrt(sum(t * t))
      t <- t / tnorm
      r <- r / tnorm

      p <- crossprod(X, t)
      c <- crossprod(Y, t)
      u <- Y %*% c
      v <- p

      if (a > 1) {
         v <- v - V %*% crossprod(V, p)
         u <- u - TT %*% crossprod(TT, u)
      }

      v <- v / sqrt(sum(v * v))

      R[, a] <- r
      V[, a] <- v
      P[, a] <- p
      TT[, a] <- t
      U[, a] <- u
      C[, a] <- c

      M <- M - tcrossprod(p)
      S <- S - v %*% crossprod(v, S)
   }

   return(list(R = R, P = P, C = C))
}

