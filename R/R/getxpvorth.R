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
