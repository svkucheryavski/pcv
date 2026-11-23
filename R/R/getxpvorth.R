#' Compute sign-flip diagonal matrix for direction correction of local loadings
#'
#' @description
#' Compares column directions of two loading matrices and returns a diagonal matrix
#' of +1/-1 signs. If the angle between corresponding columns exceeds 90 degrees,
#' the sign is -1 (flip), otherwise +1 (keep). Input dot products are clamped to
#' [-1, 1] before \code{acos} to avoid NaN from floating-point overshoot.
#'
#' @param V1
#' matrix with reference vectors as columns (e.g. global loadings)
#' @param V2
#' matrix with vectors to align as columns (e.g. local loadings)
#'
#' @return
#' Diagonal matrix (ncol x ncol) with +1 or -1 on the diagonal
#'
#' @keywords internal
flipSigns <- function(V1, V2) {
   cosines <- pmin(pmax(colSums(V1 * V2), -1), 1)
   a <- acos(cosines) < pi / 2
   diag(a * 2 - 1, ncol(V1), ncol(V1))
}


#' Generates the orthogonal part for Xpv
#' @param q.k
#' vector with orthogonal distances for cross-validation set for segment k
#' @param X.k
#' matrix with local validation set for segment k
#' @param PRM
#' projecton matrix for orthogonalization of residuals
#'
#' @return
#' A matrix with orthogonal part for Xpv
#'
#' @importFrom stats runif
#'
#' @keywords internal
getxpvorth <- function(q.k, X.k, PRM) {

   nobj <- length(q.k)

   # compute the orthogonal component of Xpv
   Xpv.orth <- matrix(runif(nobj^2, -1, 1), nobj, nobj) %*% X.k                # - project Xk to a random vector
   n1 <- sqrt(colSums(Xpv.orth^2))                                             # - compute Eucledian norm for columns
   n1[n1 < (100 *.Machine$double.eps)] <- 1.0                                  # - to avoid division by zero
   Xpv.orth <- Xpv.orth %*% diag(1 / n1)                                       # - normalize columns
   Xpv.orth <- Xpv.orth %*% PRM                                                # - orthogonalize to global component space
   n2 <- rowSums(Xpv.orth^2)                                                   # - compute sum of squared values for each row
   n2[n2 < (100 *.Machine$double.eps)] <- 1.0                                  # - to avoid division by zero
   Xpv.orth <- diag(sqrt(q.k / n2), nobj, nobj) %*% Xpv.orth  # - rescale rows

   return (Xpv.orth)
}
