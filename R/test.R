# Simple tests for the PCV R implementation
rm(list = ls())
source("pcv.R")
I <- 100
J <- 50
A <- 10
K <- 4

X <- matrix(rnorm(I * J), I, J)

params <- list()
params[[1]] <- list(X = X)
params[[2]] <- list(X = X, ncomp = 1)
params[[3]] <- list(X = X, ncomp = A)
params[[4]] <- list(X = X, ncomp = A, nseg = 4)
params[[5]] <- list(X = X, ncomp = A, nseg = 10)
params[[6]] <- list(X = X, ncomp = A, nseg = nrow(X))
params[[7]] <- list(X = X, ncomp = A, nseg = 4, scale = TRUE)
params[[8]] <- list(X = X, ncomp = A, nseg = 10, scale = TRUE)
params[[9]] <- list(X = X, ncomp = A, nseg = nrow(X), scale = TRUE)

cat("\nRunning: ")
for (i in seq_along(params)) {
   cat(i, " ")
   Xpv <- do.call(pcv, params[[i]])
   stopifnot(nrow(Xpv) == nrow(X))
   stopifnot(ncol(Xpv) == ncol(X))
   stopifnot(ks.test(X, Xpv)$p.value > 0.01)
}
cat("\n")

