#' Procrustes cross-validation for PLS models
#'
#' @param X
#' matrix with predictors from the training set.
#' @param y
#' vector with response values from the training set.
#' @param ncomp
#' number of components to use (more than the expected optimal number).
#' @param center
#' logical, to center or not the data sets
#' @param scale
#' logical, to scale or not the data sets
#' @param cv
#' which split method to use for cross-validation (see description for details).
#'
#' @description
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
#' Pseudo-validation matrix (same size as X)
#'
#' @export
pcvpls <- function(X, y, ncomp = min(nrow(X) - 1, col(X), 30), center = TRUE, scale = FALSE, cv = list("ven", 4)) {

   funlist <- list(

      # computes global PLS model
      getglobalmodel = function(X, y, ncomp) {

         m <- simpls(X, y, ncomp)
         P <- m$P
         C <- m$C
         R <- m$R
         Pi <- (eye(ncol(X)) - tcrossprod(R, P))
         return(list(P = P, Pi = Pi, C = C, R = R, ncomp = ncomp))
      },

      # computes local PLS model
      getlocalmodel = function(X.c, y.c, m) {

         # local model
         m.k <- simpls(X.c, y.c, m$ncomp)
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

   return(pcvreg(X, y, ncomp, cv = cv, center = center, scale = scale, funlist = funlist))
}


#' Procrustes cross-validation for PCR models
#'
#' @param X
#' matrix with predictors from the training set.
#' @param y
#' vector with response values from the training set.
#' @param ncomp
#' number of components to use (more than the expected optimal number).
#' @param center
#' logical, to center or not the data sets
#' @param scale
#' logical, to scale or not the data sets
#' @param cv
#' which split method to use for cross-validation (see description of method `pcvpls()` for details).
#'
#' @description
#' The method computes pseudo-validation matrix Xpv, based on PCR decomposition of calibration
#' set {X, y} and cross-validation.
#'
#' @return
#' Pseudo-validation matrix (same size as X)
#'
#' @export
pcvpcr <- function(X, y, ncomp = min(nrow(X) - 1, col(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE) {

   funlist <- list(

      # computes global PCR model
      getglobalmodel = function(X, y, ncomp) {
         P <- svd(X)$v[, 1:ncomp, drop = FALSE]
         T <- X %*% P
         C <- t(solve(crossprod(T)) %*% crossprod(T, y))
         I <- eye(ncol(X))
         Pi <- I - tcrossprod(P)
         return(list(P = P, C = C, Pi = Pi, ncomp = ncomp))
      },

      # computes local PCR model
      getlocalmodel = function(X.c, y.c, m) {
         P.k <- svd(X.c)$v[, seq_len(ncomp), drop = FALSE]
         aa <- acos(colSums(m$P * P.k)) < (pi / 2)
         P.k <- P.k %*% diag(aa * 2 - 1, ncol(P.k), ncol(P.k))

         T.c <- X.c %*% P.k
         C.k <- t(solve(crossprod(T.c)) %*% crossprod(T.c, y.c))

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

   return(pcvreg(X, y, ncomp, cv = cv, center = center, scale = scale, funlist = funlist))
}


#' Procrustes cross-validation for multivariate regression models
#'
#' @description
#' This is a generic method, use `pcvpls()` or `pcvpcr()` instead.
#'
pcvreg <- function(X, Y, ncomp = min(nrow(X) - 1, col(X), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE, funlist = list()) {

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
   ind <- crossval(cv, nrow(X), Y[, 1])

   # get number of segments
   nseg <- max(ind)

   # correct maximum number of components
   max.ncomp <- min(nrow(X) - ceiling(nrow(X) / nseg) - 1 , ncol(X))
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
      print(which(ind.k))
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


pcvpca <- function(X, ncomp = min(nrow(X) - 1, col(X), 20), cv = list("ven", 4),
   center = TRUE, scale = FALSE) {

   # keep names if any
   attrs <- attributes(X)

   mX <- apply(X, 2, mean)
   sX <- if (scale) apply(X, 2, sd) else rep(1, ncol(X))

   # autoscale the calibration set
   X <- scale(X, center = mX, scale = sX)

   # indices for cross-validation
   ind <- crossval(cv, nrow(X), seq_len(nrow(X)))

   # get number of segments
   nseg <- max(ind)

   # correct maximum number of components
   max.ncomp <- min(nrow(X) - ceiling(nrow(X) / nseg) - 1 , ncol(X))
   if (ncomp > max.ncomp) {
      ncomp <- max.ncomp
   }

   # create a global model
   P <- svd(X)$v[, seq_len(ncomp), drop = FALSE]
   Pi <- eye(nrow(P)) - tcrossprod(P)

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

################################################
# Helper functions                             #
################################################

#' Generates the orthogonal part for Xpv
#' @param q.k
#' vector with orthogonal distances for cross-validation set
#' @param PRM
#' projecton matrix for orthogonalization of residuals
getxpvorth <- function(q.k, X.k, PRM) {

   nobj <- length(q.k)

   # compute the orthogonal component of Xpv
   Xpv.orth <- matrix(runif(nobj^2, -1, 1), nobj, nobj) %*% X.k           # - project Xk to a random vector
   Xpv.orth <- Xpv.orth %*% diag(1 / sqrt(colSums(Xpv.orth^2)))                # - normalize columns
   Xpv.orth <- Xpv.orth %*% PRM                                                # - orthogonalize to global component space
   Xpv.orth <- diag(sqrt(q.k / rowSums(Xpv.orth^2)), nobj, nobj) %*% Xpv.orth  # - rescale rows

   return (Xpv.orth)
}

#' Generate sequence of indices for cross-validation
#'
#' @description
#' Generates and returns sequence of object indices for each segment in random segmented
#' cross-validation
#'
#' @param nobj
#' number of objects in a dataset
#' @param cv
#' cross-validation settings, can be a number or a list. If cv is a number, it will be
#' used as a number of segments for random cross-validation (if cv = 1, full cross-validation
#' will be preformed), if it is a list, the following syntax can be used:
#' cv = list('rand', nseg, nrep) for random repeated cross-validation with nseg segments and nrep
#' repetitions or cv = list('ven', nseg) for systematic splits to nseg segments ('venetian blinds').
#' @param resp
#' vector with response values to use in case of venetian blinds
#'
#' @return
#' matrix with object indices for each segment
#'
#' @export
crossval <- function(cv = 1, nobj = NULL, resp = NULL) {

   # get cross-validation parameters
   if (is.null(nobj)) nobj <- length(resp)

   # if user already provided matrix with values - return it
   if (is.numeric(cv) && length(cv) == nobj) return(as.matrix(cv))

   p <- getCrossvalParams(cv = cv, nobj = nobj)
   if (!(p$type %in% c("rand", "ven", "loo"))) {
      stop("Wrong name for cross-validation method.")
   }

   # check number of repetitions
   if (p$nrep < 1 || p$nrep > 100) {
      stop("Wrong value for cv repetitions (should be between 1 and 100).")
   }

   # check number of segments
   if (p$nseg < 2 || p$nseg > nobj) {
      stop("Wrong value for number of segments (should be between 2 and number of objects).")
   }

   if (p$type == "loo") {
      return(matrix(seq_len(nobj), ncol = 1))
   }

   if (p$type == "rand") {
      return(sapply(seq_len(p$nrep), function(i) rep(seq_len(p$nseg), length.out = nobj)[sample(nobj)]))
   }

   if (p$type == "ven") {
      ind <- if (is.null(resp)) seq_len(nobj) else order(resp)
      return(matrix(rep(seq_len(p$nseg), length.out = nobj)[ind], ncol = 1))
   }

   stop("Something went wrong.")
}


#' Returns parameters for cross-validation based on 'cv' value
#'
#' @param cv
#' settings for cross-validation provided by user
#' @param nobj
#' number of objects in calibration set
#'
getCrossvalParams <- function(cv, nobj) {

   nrep <- 1

   # random
   if (is.numeric(cv)) {
      return(
         list(
            type = "rand",
            nrep = 1,
            nseg = if (cv == 1) nobj else cv
         )
      )
   }

   # leave one out
   type <- cv[[1]]
   if (type == "loo") {
      return(
         list(
            type = "loo",
            nrep = 1,
            nseg = nobj
         )
      )
   }

   # venetian blinds
   nseg <- cv[[2]]
   if (type == "ven") {
      return(
         list(
            type = "ven",
            nrep = nrep,
            nseg = nseg
         )
      )
   }

   nrep <- if (length(cv) == 3) cv[[3]] else 1
   return(
      list(
         type = type,
         nrep = nrep,
         nseg = nseg
      )
   )
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
#' @param x
#' a matrix with x values (predictors)
#' @param y
#' a matrix with y values (responses)
#' @param ncomp
#' number of components to calculate
#' @param cv
#' logical, is model calibrated during cross-validation or not
#'
#' @return
#' a list with computed regression coefficients, loadings and scores for x and y matrices,
#' and weights.
#'
#' @references
#' [1]. S. de Jong. SIMPLS: An Alternative approach to partial least squares regression.
#' Chemometrics and Intelligent Laboratory Systems, 18, 1993 (251-263).
#'
simpls <- function(x, y, ncomp) {

   x <- as.matrix(x)
   y <- as.matrix(y)

   nobj  <- nrow(x)
   npred <- ncol(x)
   nresp <- ncol(y)

   # initial estimation
   S <- crossprod(x, y)
   M <- crossprod(x)

   # prepare space for results
   C <- matrix(0, nrow = nresp, ncol = ncomp)
   R <- V <- P <- matrix(0, nrow = npred, ncol = ncomp)
   TT <- U <- matrix(0, nrow = nobj, ncol = ncomp)


   # loop for each components
   for (a in seq_len(ncomp)) {

      r <- svd(S, nu = 1, nv = 0)$u
      t <- x %*% r

      tnorm <- sqrt(sum(t * t))
      t <- t / tnorm
      r <- r / tnorm

      p <- crossprod(x, t)
      c <- crossprod(y, t)
      u <- y %*% c
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

#' Create the identity matrix
#'
#' @param n
#' Size of the matrix
#'
#' @return
#' The identity matrix (n x n)
#'
#' @export
eye <- function(n) {
   X <- matrix(0, n, n)
   diag(X) <- 1
   return(X)
}
