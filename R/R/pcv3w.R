#' Procrustes cross-validation for three-way data: full unfolding (FU)
#'
#' @param X3w
#' three-way array with the calibration set, dimension \code{c(I, J, K)} (samples x mode 2 x mode 3).
#' @param ncomp
#' number of PCA components (more than the expected optimal number of components of the
#' three-way model, see details).
#' @param cv
#' which split method to use for cross-validation (see \code{\link{pcvpca}} for details).
#' @param center
#' logical, to center the unfolded data or not.
#' @param scale
#' logical, to scale the unfolded data or not (each of the J x K columns separately).
#' @param seed
#' optional integer, makes random splits and random orthogonal parts reproducible. The random
#' number generator state of the caller is restored afterwards.
#'
#' @details
#' The method unfolds the three-way array to a matrix with I rows and J*K columns, generates
#' PV-set for the matrix using \code{\link{pcvpca}} with global CV scope, and folds the result
#' back to a three-way array.
#'
#' The PV-set reproduces cross-validation of the PCA model of the unfolded data. Use
#' enough components, so the PCA subspace contains the subspace of the three-way models
#' (e.g. PARAFAC or Tucker3) you are going to validate.
#'
#' @return
#' Three-way array with PV-set, same dimension as \code{X3w}.
#'
#' @seealso
#' \code{\link{pcv3wbu}}, \code{\link{pcv3whsv}}
#'
#' @export
pcv3wfu <- function(X3w, ncomp = min(dim(X3w)[1] - 1, prod(dim(X3w)[2:3]), 30), cv = list("ven", 4),
   center = TRUE, scale = FALSE, seed = NULL) {

   d <- getdim3w(X3w)

   # unfold to I x (J*K), mode 2 runs fastest, generate PV-set and fold back
   X <- matrix(X3w, d[1], d[2] * d[3])
   withseed3w(seed, {
      cvind <- splitcv3w(cv, d[1], "pcv3wfu")
      Xpv <- pcvpca(X, ncomp = ncomp, cv = cvind, center = center, scale = scale, cv.scope = "global")
   })

   Xpv3w <- array(Xpv, d)
   dimnames(Xpv3w) <- dimnames(X3w)
   return(Xpv3w)
}


#' Procrustes cross-validation for three-way data: block unfolding (BU)
#'
#' @param X3w
#' three-way array with the calibration set, dimension \code{c(I, J, K)} (samples x mode 2 x mode 3).
#' @param mode
#' which mode to split into blocks (2 or 3).
#' @param ncomp
#' number of PCA components for each block.
#' @param cv
#' which split method to use for cross-validation (see \code{\link{pcvpca}} for details).
#' @param center
#' logical, to center each block or not.
#' @param scale
#' logical, to scale each block or not.
#' @param seed
#' optional integer, makes random splits and random orthogonal parts reproducible. The random
#' number generator state of the caller is restored afterwards.
#'
#' @details
#' The method splits the three-way array into J matrices of size I x K (if \code{mode = 2}) or
#' K matrices of size I x J (if \code{mode = 3}). Then it generates PV-set for each matrix using
#' \code{\link{pcvpca}} with global CV scope and combines the PV-sets into a three-way array.
#' All blocks use the same split of the samples into segments.
#'
#' @return
#' Three-way array with PV-set, same dimension as \code{X3w}. Attribute \code{cvind} contains
#' the segment index of each sample.
#'
#' @seealso
#' \code{\link{pcv3wfu}}, \code{\link{pcv3whsv}}
#'
#' @export
pcv3wbu <- function(X3w, mode = 2, ncomp = min(dim(X3w)[1] - 1, dim(X3w)[if (mode == 2) 3 else 2], 30),
   cv = list("ven", 4), center = TRUE, scale = FALSE, seed = NULL) {

   d <- getdim3w(X3w)
   if (length(mode) != 1 || !(mode %in% c(2, 3))) stop("pcv3wbu: 'mode' must be 2 or 3.")

   withseed3w(seed, {
      cvind <- splitcv3w(cv, d[1], "pcv3wbu")
      Xpv3w <- array(0, dim = d)
      for (b in seq_len(d[mode])) {
         X <- if (mode == 2) X3w[, b, ] else X3w[, , b]
         Xpv <- pcvpca(X, ncomp = ncomp, cv = cvind, center = center, scale = scale,
            cv.scope = "global")
         if (mode == 2) Xpv3w[, b, ] <- Xpv else Xpv3w[, , b] <- Xpv
      }
   })

   dimnames(Xpv3w) <- dimnames(X3w)
   attr(Xpv3w, "cvind") <- cvind
   return(Xpv3w)
}


#' Procrustes cross-validation for three-way data: HO-SVD based (HSV)
#'
#' @param X3w
#' three-way array with the calibration set, dimension \code{c(I, J, K)} (samples x mode 2 x mode 3),
#' all values must be finite.
#' @param ncomp
#' vector with two positive integers, \code{c(R2, R3)}, number of components in modes 2 and 3.
#' @param cv
#' which split method to use for cross-validation (see \code{\link{pcvpca}} for details).
#' @param center
#' logical, to subtract the mean slab or not.
#' @param scale
#' scaling within modes: \code{"none"} (or \code{FALSE}), \code{"mode2"}, \code{"mode3"} or
#' \code{"both"} (or \code{TRUE}), see details. A variable with zero variance (also when it is
#' constant over the samples) cannot be scaled and gives an error.
#' @param orth
#' how to construct the orthogonal part: \code{"random"} (random mix of the segment's slabs,
#' similar to the two-way method) or \code{"canonical"} (each slab's own local residual, no
#' randomness).
#' @param align
#' how to match local components to the global ones: \code{"optimal"} (maximal sum of absolute
#' cosines) or \code{"greedy"} (global components in order of variance choose first).
#' @param seed
#' optional integer, makes random splits and random orthogonal parts reproducible. The random
#' number generator state of the caller is restored afterwards.
#'
#' @details
#' The method models the three-way array with truncated HO-SVD in modes 2 and 3 (a Tucker2 model)
#' and generates PV-set, which, when projected on the global model, gives the same scores and
#' orthogonal distances as cross-validation of the model, sample by sample.
#'
#' Scaling divides mode-2 and/or mode-3 variables by their standard deviation after centering.
#' With \code{"both"}, mode 2 is scaled first and mode 3 is then scaled using the mode 2 scaled
#' data. The PV-set is returned in the original units.
#'
#' A PV-set generated with \code{c(R2, R3)} components also reproduces cross-validation of the
#' models with fewer components, \code{r2 <= R2} and \code{r3 <= R3}, unless the local components
#' are matched to the global ones across the boundary of the first \code{r2} (or \code{r3})
#' components. The numbers of components for which this holds are returned in attribute
#' \code{valid.subranks}, and the function gives a warning if some of them are missing.
#'
#' @return
#' Three-way array with PV-set, same dimension as \code{X3w}, with attributes:
#' \code{ncomp} (number of components actually used), \code{cvind} (segment index of each
#' sample), \code{orth}, \code{align}, \code{scale}, \code{prep} (global preprocessing
#' parameters), \code{alignment} (matched local component and its absolute cosine for every
#' segment, mode and global component) and \code{valid.subranks}.
#'
#' @seealso
#' \code{\link{pcv3wfu}}, \code{\link{pcv3wbu}}
#'
#' @export
pcv3whsv <- function(X3w, ncomp, cv = list("ven", 4), center = TRUE, scale = FALSE,
   orth = c("random", "canonical"), align = c("optimal", "greedy"), seed = NULL) {

   orth <- match.arg(orth)
   align <- match.arg(align)
   scale <- normalizescale3w(scale)
   d <- getdim3w(X3w)

   if (!all(is.finite(X3w))) {
      stop("pcv3whsv: 'X3w' contains NA/NaN/Inf values, interpolate or impute them first.")
   }

   if (length(ncomp) != 2 || !is.numeric(ncomp) || anyNA(ncomp) || any(ncomp < 1) ||
         any(ncomp != round(ncomp))) {
      stop("pcv3whsv: 'ncomp' must be two positive integers, c(R2, R3).")
   }

   withseed3w(seed, {
      cvind <- splitcv3w(cv, d[1], "pcv3whsv")
      res <- pcv3whsv.run(X3w, ncomp, cvind, center, scale, orth, align)
   })

   return(res)
}


#' Main loop of pcv3whsv (arguments already checked)
#' @keywords internal
pcv3whsv.run <- function(X3w, ncomp, cvind, center, scale, orth, align) {

   d <- dim(X3w)
   nseg <- max(cvind)

   # global preprocessing and global factors
   pg <- prep3w(X3w, center, scale)
   gmf <- modefactors3w(pg$X3wp, ncomp[1], ncomp[2])
   r2min <- gmf$r2
   r3min <- gmf$r3

   # local factors for every segment
   lmf <- vector("list", nseg)
   for (k in seq_len(nseg)) {
      lmf[[k]] <- modefactors3w(pg$X3wp[cvind != k, , , drop = FALSE], ncomp[1], ncomp[2])
      r2min <- min(r2min, lmf[[k]]$r2)
      r3min <- min(r3min, lmf[[k]]$r3)
   }

   # ranks attainable by the global and all local models
   R2 <- min(ncomp[1], r2min)
   R3 <- min(ncomp[2], r3min)

   if (R2 < 1 || R3 < 1) {
      stop("pcv3whsv: the data cannot supply at least one component in every mode for all segments.")
   }

   if (R2 < ncomp[1] || R3 < ncomp[2]) {
      warning(sprintf("pcv3whsv: 'ncomp' reduced to c(%d, %d) (maximum for all segments).", R2, R3))
   }

   U2 <- gmf$U2[, seq_len(R2), drop = FALSE]
   U3 <- gmf$U3[, seq_len(R3), drop = FALSE]
   cdir <- complementdir3w(U2, U3)

   X3wpv <- array(0, dim = d)
   adiag <- vector("list", nseg)

   for (k in seq_len(nseg)) {
      indk <- which(cvind == k)
      Xk <- lapply(indk, function(i) slab3w(pg$X3wp, i))

      # local factors matched to the global ones
      U2k <- align3w(U2, lmf[[k]]$U2[, seq_len(R2), drop = FALSE], align)
      U3k <- align3w(U3, lmf[[k]]$U3[, seq_len(R3), drop = FALSE], align)
      adiag[[k]] <- rbind(aligndiag3w(k, 2L, U2k), aligndiag3w(k, 3L, U3k))

      # explained part: local scores with global factors; residuals of the local model
      Xhat <- vector("list", length(Xk))
      Eres <- vector("list", length(Xk))
      for (m in seq_along(Xk)) {
         Ti <- slabscores3w(Xk[[m]], U2k, U3k)
         Eres[[m]] <- Xk[[m]] - U2k %*% Ti %*% t(U3k)
         Xhat[[m]] <- U2 %*% Ti %*% t(U3)
      }

      Eo <- getxpvorth3w(Xk, Eres, U2, U3, cdir, orth)
      for (m in seq_along(indk)) {
         X3wpv[indk[m], , ] <- Xhat[[m]] + Eo[[m]]
      }
   }

   X3wpv <- prep3w.remove(X3wpv, pg)
   dimnames(X3wpv) <- dimnames(X3w)

   # sub-rank r is valid if in every segment components 1..r are matched within 1..r
   adiag <- do.call(rbind, adiag)
   rownames(adiag) <- NULL
   validsub <- function(mode, R) {
      ok <- rep(TRUE, R)
      for (k in unique(adiag$segment)) {
         perm <- adiag$local.col[adiag$segment == k & adiag$mode == mode]
         ok <- ok & (cummax(perm) == seq_len(R))
      }
      which(ok)
   }
   v2 <- validsub(2L, R2)
   v3 <- validsub(3L, R3)

   if (length(v2) < R2 || length(v3) < R3) {
      warning("pcv3whsv: PV-set does not reproduce cross-validation for some smaller numbers of ",
         "components, see attribute 'valid.subranks'.")
   }

   attr(X3wpv, "ncomp") <- c(R2, R3)
   attr(X3wpv, "cvind") <- cvind
   attr(X3wpv, "orth") <- orth
   attr(X3wpv, "align") <- align
   attr(X3wpv, "scale") <- scale
   attr(X3wpv, "prep") <- pg[c("M", "S2", "S3", "center", "scale")]
   attr(X3wpv, "alignment") <- adiag
   attr(X3wpv, "valid.subranks") <- list(r2 = v2, r3 = v3)

   return(X3wpv)
}


###################################################
# Checks and settings                             #
###################################################

#' Check that X3w is a numeric three-way array and return its dimension
#' @keywords internal
getdim3w <- function(X3w) {
   d <- dim(X3w)
   if (!is.numeric(X3w) || length(d) != 3) stop("'X3w' must be a numeric three-way array.")
   if (any(d < 2)) stop("'X3w' must have at least 2 values in each mode.")
   return(d)
}


#' Convert the 'scale' argument to one of "none", "mode2", "mode3", "both"
#' @keywords internal
normalizescale3w <- function(scale) {
   if (isFALSE(scale)) return("none")
   if (isTRUE(scale)) return("both")
   if (!is.character(scale) || length(scale) != 1) {
      stop("'scale' must be FALSE, TRUE or one of \"none\", \"mode2\", \"mode3\", \"both\".")
   }
   match.arg(scale, c("none", "mode2", "mode3", "both"))
}


#' Segment number of each sample (first repetition only), checked to be 1, 2, ..., nseg
#' @keywords internal
splitcv3w <- function(cv, nobj, fname) {
   cvfull <- pcvcrossval(cv, nobj)
   if (ncol(cvfull) > 1) {
      warning(fname, ": 'cv' specifies several repetitions, only the first one is used.")
   }

   cvind <- cvfull[, 1]
   nseg <- if (is.numeric(cvind) && !anyNA(cvind)) max(cvind) else NA
   if (is.na(nseg) || any(cvind != round(cvind)) || nseg < 2 || any(cvind < 1) ||
         !all(seq_len(nseg) %in% cvind)) {
      stop(fname, ": segment numbers in 'cv' must be 1, 2, ..., n (n >= 2), each used at least once.")
   }

   return(as.integer(cvind))
}


#' Evaluate code with a given seed and restore the caller's random number generator state
#' @keywords internal
withseed3w <- function(seed, code) {
   if (is.null(seed)) return(code)

   genv <- globalenv()
   hasseed <- exists(".Random.seed", envir = genv, inherits = FALSE)
   oldseed <- if (hasseed) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
   on.exit(
      if (hasseed) assign(".Random.seed", oldseed, envir = genv)
      else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
   )

   set.seed(seed)
   return(code)
}


###################################################
# Preprocessing                                   #
###################################################

#' Center (mean slab) and scale (within modes) a three-way array
#'
#' @return list with preprocessed array \code{X3wp} and parameters \code{M, S2, S3, center, scale}
#' @keywords internal
prep3w <- function(X3w, center = TRUE, scale = "none") {

   d <- dim(X3w)
   M <- matrix(0, d[2], d[3])
   S2 <- rep(1, d[2])
   S3 <- rep(1, d[3])

   # zero variance check: sd below the rounding noise of centering (10 * I * eps times the
   # largest absolute value of the variable before centering) counts as zero
   A <- abs(X3w)
   tol <- 10 * d[1] * .Machine$double.eps

   X3wp <- X3w
   if (center) {
      M <- apply(X3wp, c(2, 3), mean)
      X3wp <- sweep(X3wp, c(2, 3), M, "-")
   }

   scale <- normalizescale3w(scale)

   if (scale %in% c("mode2", "both")) {
      S2 <- apply(X3wp, 2, sd)
      bad <- which(!is.finite(S2) | !(S2 > tol * apply(A, 2, max)))
      if (length(bad) > 0) {
         stop("prep3w: cannot scale mode 2, zero variance in variable(s) ",
            paste(bad[seq_len(min(5, length(bad)))], collapse = ", "), ".")
      }
      X3wp <- sweep(X3wp, 2, S2, "/")
      A <- sweep(A, 2, S2, "/")
   }

   if (scale %in% c("mode3", "both")) {
      S3 <- apply(X3wp, 3, sd)
      bad <- which(!is.finite(S3) | !(S3 > tol * apply(A, 3, max)))
      if (length(bad) > 0) {
         stop("prep3w: cannot scale mode 3, zero variance in variable(s) ",
            paste(bad[seq_len(min(5, length(bad)))], collapse = ", "), ".")
      }
      X3wp <- sweep(X3wp, 3, S3, "/")
   }

   return(list(X3wp = X3wp, M = M, S2 = S2, S3 = S3, center = center, scale = scale))
}


#' Undo preprocessing made by prep3w (back to original units)
#' @keywords internal
prep3w.remove <- function(X3wp, p) {
   X3w <- X3wp
   if (p$scale %in% c("mode3", "both")) X3w <- sweep(X3w, 3, p$S3, "*")
   if (p$scale %in% c("mode2", "both")) X3w <- sweep(X3w, 2, p$S2, "*")
   if (p$center) X3w <- sweep(X3w, c(2, 3), p$M, "+")
   return(X3w)
}


###################################################
# Model helpers                                   #
###################################################

#' Number of singular values above the relative tolerance
#' @keywords internal
numrank3w <- function(d, tol = sqrt(.Machine$double.eps)) {
   if (length(d) == 0 || d[1] <= 0) return(0L)
   sum(d > tol * d[1])
}


#' Slab i of a three-way array as a J x K matrix
#' @keywords internal
slab3w <- function(X3w, i) {
   d <- dim(X3w)
   matrix(X3w[i, , , drop = FALSE], d[2], d[3])
}


#' Mode-n unfolding (mode 2: J x (K*I), mode 3: K x (J*I))
#' @keywords internal
unfold3w <- function(X3w, mode) {
   d <- dim(X3w)
   if (mode == 2) return(matrix(aperm(X3w, c(2, 3, 1)), nrow = d[2]))
   if (mode == 3) return(matrix(aperm(X3w, c(3, 2, 1)), nrow = d[3]))
   stop("unfold3w: mode must be 2 or 3.")
}


#' Mode 2 and mode 3 factors (truncated SVD of unfoldings) and numerical ranks
#'
#' @return list with \code{U2, U3, r2, r3}
#' @keywords internal
modefactors3w <- function(X3w, R2, R3) {
   M2 <- unfold3w(X3w, 2)
   M3 <- unfold3w(X3w, 3)
   s2 <- svd(M2, nu = min(R2, nrow(M2), ncol(M2)), nv = 0)
   s3 <- svd(M3, nu = min(R3, nrow(M3), ncol(M3)), nv = 0)
   list(U2 = s2$u, U3 = s3$u, r2 = numrank3w(s2$d), r3 = numrank3w(s3$d))
}


#' Two-sided scores of a slab: t(U2) X U3
#' @keywords internal
slabscores3w <- function(Xi, U2, U3) crossprod(U2, Xi) %*% U3


#' Projection of a slab onto range(U2) x range(U3)
#' @keywords internal
proj3w <- function(S, U2, U3) U2 %*% slabscores3w(S, U2, U3) %*% t(U3)


###################################################
# Orthogonal part                                 #
###################################################

#' Unit vector orthogonal to the columns of U
#' @keywords internal
complvec3w <- function(U) {
   j <- which.min(rowSums(U^2))
   e <- numeric(nrow(U))
   e[j] <- 1
   w <- e - U %*% crossprod(U, e)
   w / sqrt(sum(w^2))
}


#' Fixed unit direction outside the model space (NULL if there is none)
#' @keywords internal
complementdir3w <- function(U2, U3) {
   # first factor with sign fixed, so the direction does not depend on the SVD sign convention
   first <- function(U) U[, 1, drop = FALSE] * sign(U[which.max(abs(U[, 1])), 1])
   if (ncol(U2) < nrow(U2)) return(complvec3w(U2) %*% t(first(U3)))
   if (ncol(U3) < nrow(U3)) return(first(U2) %*% t(complvec3w(U3)))
   return(NULL)
}


#' Orthogonal parts of the PV-set slabs of one segment
#'
#' @description
#' Each orthogonal part lies outside the global model space and has the same norm as the
#' local residual of the slab. Its direction comes from a random mix of the segment's slabs
#' (\code{"random"}) or from the slab's own local residual (\code{"canonical"}).
#'
#' @param Xk list of held-out slabs of the segment (preprocessed)
#' @param Eres list of local residuals of the slabs
#' @param U2,U3 global factors
#' @param cdir fallback direction from \code{complementdir3w}
#' @param orth "random" or "canonical"
#'
#' @return list of J x K matrices
#' @keywords internal
getxpvorth3w <- function(Xk, Eres, U2, U3, cdir, orth = c("random", "canonical")) {
   orth <- match.arg(orth)
   nk <- length(Xk)
   J <- nrow(U2)
   K <- nrow(U3)
   eps <- 100 * .Machine$double.eps

   Z <- if (orth == "random") matrix(stats::runif(nk * nk, -1, 1), nk, nk) else NULL
   E <- vector("list", nk)

   for (i in seq_len(nk)) {
      qi <- sum(Eres[[i]]^2)
      if (qi <= 0) {
         E[[i]] <- matrix(0, J, K)
         next
      }

      if (orth == "canonical") {
         src <- Eres[[i]]
      } else {
         src <- matrix(0, J, K)
         for (j in seq_len(nk)) src <- src + Z[i, j] * Xk[[j]]
      }

      Eo <- src - proj3w(src, U2, U3)

      # source lies (almost) inside the model space: use the fixed direction
      if (sum(Eo^2) <= eps * sum(src^2)) {
         if (is.null(cdir)) {
            E[[i]] <- matrix(0, J, K)
            next
         }
         Eo <- cdir
      }

      E[[i]] <- Eo * sqrt(qi / sum(Eo^2))
   }

   return(E)
}


###################################################
# Matching local and global components            #
###################################################

#' Assignment maximizing the sum of M[i, p[i]] (Hungarian method)
#'
#' @return integer vector p, p[i] is the column assigned to row i
#' @keywords internal
lsap3w <- function(M) {
   n <- nrow(M)
   C <- max(M) - M

   # potentials and matching use index 0 as a sentinel, stored at position 1
   u <- numeric(n + 1)
   v <- numeric(n + 1)
   p <- integer(n + 1)
   way <- integer(n + 1)

   for (i in seq_len(n)) {
      p[1] <- i
      j0 <- 0L
      minv <- rep(Inf, n + 1)
      used <- rep(FALSE, n + 1)
      repeat {
         used[j0 + 1] <- TRUE
         i0 <- p[j0 + 1]
         delta <- Inf
         j1 <- 0L
         for (j in seq_len(n)) {
            if (used[j + 1]) next
            cur <- C[i0, j] - u[i0 + 1] - v[j + 1]
            if (cur < minv[j + 1]) {
               minv[j + 1] <- cur
               way[j + 1] <- j0
            }
            if (minv[j + 1] < delta) {
               delta <- minv[j + 1]
               j1 <- j
            }
         }
         for (j in 0:n) {
            if (used[j + 1]) {
               u[p[j + 1] + 1] <- u[p[j + 1] + 1] + delta
               v[j + 1] <- v[j + 1] - delta
            } else {
               minv[j + 1] <- minv[j + 1] - delta
            }
         }
         j0 <- j1
         if (p[j0 + 1] == 0) break
      }
      repeat {
         j1 <- way[j0 + 1]
         p[j0 + 1] <- p[j1 + 1]
         j0 <- j1
         if (j0 == 0) break
      }
   }

   res <- integer(n)
   for (j in seq_len(n)) res[p[j + 1]] <- j
   return(res)
}


#' Match local factors to global ones by permutation and sign
#'
#' @param Ug global factors (columns are components)
#' @param Ul local factors (same number of columns)
#' @param method "optimal" or "greedy"
#'
#' @return matched local factors with attributes \code{perm} (local column for each global
#' component) and \code{cos} (absolute cosines of the matched pairs)
#' @keywords internal
align3w <- function(Ug, Ul, method = c("optimal", "greedy")) {
   method <- match.arg(method)
   R <- ncol(Ug)
   if (ncol(Ul) != R) stop("align3w: factor matrices must have the same number of columns.")

   M <- abs(crossprod(Ug, Ul))
   if (method == "optimal") {
      assigned <- lsap3w(M)
   } else {
      assigned <- integer(R)
      usedlocal <- logical(R)
      for (g in seq_len(R)) {
         cand <- which(!usedlocal)
         best <- cand[which.max(M[g, cand])]
         assigned[g] <- best
         usedlocal[best] <- TRUE
      }
   }

   Ual <- Ul[, assigned, drop = FALSE]
   sgn <- sign(diag(crossprod(Ug, Ual)))
   sgn[sgn == 0] <- 1

   structure(Ual %*% diag(sgn, R, R), perm = assigned, cos = M[cbind(seq_len(R), assigned)])
}


#' Matching diagnostics of one segment and mode as a data frame
#' @keywords internal
aligndiag3w <- function(segment, mode, aligned) {
   perm <- attr(aligned, "perm")
   data.frame(segment = segment, mode = mode, component = seq_along(perm),
      local.col = perm, cosine = attr(aligned, "cos"))
}
