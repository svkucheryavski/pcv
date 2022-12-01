# Simple tests for the PCV R implementation
rm(list = ls())
source("pcv.R")
I <- 100
J <- 50
A <- 30
K <- 4

X <- matrix(rnorm(I * J), I, J)

params <- list()
params[[1]] <- list(X = X)
params[[2]] <- list(X = X, ncomp = 1)
params[[3]] <- list(X = X, ncomp = A)
params[[4]] <- list(X = X, ncomp = A, cv = 10)
params[[5]] <- list(X = X, ncomp = A, cv = 10, scale = TRUE)
params[[6]] <- list(X = X, ncomp = A, cv = list("ven", 4))
params[[7]] <- list(X = X, ncomp = A, cv = list("ven", 4), scale = TRUE)
params[[8]] <- list(X = X, ncomp = A, cv = list("loo"), scale = TRUE)
params[[9]] <- list(X = X, ncomp = A, cv = list("loo"))

cat("\nRunning tests for pcvpca(): ")
for (i in seq_along(params)) {
   cat(i, " ")
   Xpv <- do.call(pcvpca, params[[i]])
   stopifnot(nrow(Xpv) == nrow(X))
   stopifnot(ncol(Xpv) == ncol(X))
   stopifnot(ks.test(X, Xpv)$p.value > 0.01)
}
cat("\n")


A <- 30
y = X %*% rnorm(J)
params <- list()
params[[1]] <- list(X = X, y = y)
params[[2]] <- list(X = X, y = y, ncomp = 1)
params[[3]] <- list(X = X, y = y, ncomp = A)
params[[4]] <- list(X = X, y = y, ncomp = A, cv = 10)
params[[5]] <- list(X = X, y = y, ncomp = A, cv = 10, scale = TRUE)
params[[6]] <- list(X = X, y = y, ncomp = A, cv = list("ven", 4))
params[[7]] <- list(X = X, y = y, ncomp = A, cv = list("ven", 4), scale = TRUE)
params[[8]] <- list(X = X, y = y, ncomp = A, cv = list("loo"), scale = TRUE)
params[[9]] <- list(X = X, y = y, ncomp = A, cv = list("loo"))

cat("\nRunning tests for pcvpcr(): ")
for (i in seq_along(params)) {
   cat(i, " ")
   cv <- if (is.null(params[[i]]$cv)) list("ven", 4) else params[[i]]$cv
   Xpv <- do.call(pcvpcr, params[[i]])
   nseg <- max(crossval(cv, nobj = nrow(X)))
   stopifnot(nrow(Xpv) == nrow(X))
   stopifnot(ncol(Xpv) == ncol(X))
   stopifnot(dim(attr(Xpv, "D")) != c(A, nseg))
}
cat("\n")

cat("\nRunning tests for pcvpls(): ")
for (i in seq_along(params)) {
   cat(i, " ")
   Xpv <- do.call(pcvpls, params[[i]])
   cv <- if (is.null(params[[i]]$cv)) list("ven", 4) else params[[i]]$cv
   nseg <- max(crossval(cv, nobj = nrow(X)))
   stopifnot(nrow(Xpv) == nrow(X))
   stopifnot(ncol(Xpv) == ncol(X))
   stopifnot(dim(attr(Xpv, "D")) != c(A, nseg))
}
cat("\n")

