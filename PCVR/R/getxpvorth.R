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
#' @export
getxpvorth <- function(q.k, X.k, PRM) {

   nobj <- length(q.k)

   # compute the orthogonal component of Xpv
   Xpv.orth <- matrix(runif(nobj^2, -1, 1), nobj, nobj) %*% X.k           # - project Xk to a random vector
   Xpv.orth <- Xpv.orth %*% diag(1 / sqrt(colSums(Xpv.orth^2)))                # - normalize columns
   Xpv.orth <- Xpv.orth %*% PRM                                                # - orthogonalize to global component space
   Xpv.orth <- diag(sqrt(q.k / rowSums(Xpv.orth^2)), nobj, nobj) %*% Xpv.orth  # - rescale rows

   return (Xpv.orth)
}
