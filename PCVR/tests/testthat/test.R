####################################
# Tests for all bugs found in 2022 #
####################################

setup({
   #pdf(file = tempfile("mdatools-test-classmodel-", fileext = ".pdf"))
})

teardown({
   #dev.off()
})

context("Tests for 'pcvpca()':")

test_that("- pcvpca() works well for random data.", {
   I <- 100
   J <- 50
   A <- 30
   K <- 4

   set.seed(42)
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

   for (i in seq_along(params)) {
      expect_silent(Xpv <- do.call(pcvpca, params[[i]]))
      expect_equal(nrow(Xpv), nrow(X))
      expect_equal(ncol(Xpv), ncol(X))
      expect_true(ks.test(X, Xpv)$p.value > 0.01)
   }
})

test_that("- pcvpca() works well for Corn data.", {
   data(corn)
   X <- corn$spectra

   Xpv <- pcvpca(X, 30, center = TRUE, scale = FALSE, cv = list("ven", 4))
   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))

   mX = apply(X, 2, mean)
   Xmc = scale(X, center = mX, scale = FALSE)
   Xpvmc = scale(Xpv, center = mX, scale = FALSE)

   P = svd(Xmc)$v[, 1:20]
   T = Xmc %*% P;
   lambda = colSums(T^2 / (nrow(X) - 1))

   Tpv = Xpvmc %*% P;

   # check the distances (reproducing example from the paper)
   a = 2;
   E = Xmc - tcrossprod(T[, 1:a], P[, 1:a]);
   Epv = Xpvmc - tcrossprod(Tpv[, 1:a], P[, 1:a]);
   q = rowSums(E^2)
   qpv = rowSums(Epv^2)
   h = rowSums(T[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   hpv = rowSums(Tpv[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   expect_equal(sum(q/mean(q) > 4), 4)
   expect_equal(sum(qpv/mean(q) > 4), 4)
   expect_equal(sum(h/mean(h) > 5), 2)
   expect_equal(sum(h/mean(h) > 5), 2)

   a = 20;
   E = Xmc - tcrossprod(T[, 1:a], P[, 1:a]);
   Epv = Xpvmc - tcrossprod(Tpv[, 1:a], P[, 1:a]);
   q = rowSums(E^2)
   qpv = rowSums(Epv^2)
   h = rowSums(T[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   hpv = rowSums(Tpv[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   expect_equal(sum(q/mean(q) > 4), 0)
   expect_equal(sum(qpv/mean(q) > 4), 5)
   expect_equal(sum(h/mean(h) > 2), 3)
   expect_equal(sum(hpv/mean(h) > 2), 9)
})

context("Tests for 'pcvpcr()':")

test_that("- pcvpcr() works well for random data.", {

   I <- 100
   J <- 50
   A <- 30

   set.seed(42)
   X <- matrix(rnorm(I * J), I, J)
   Y <- X %*% runif(J)

   params <- list()
   params[[1]] <- list(X = X, Y = Y)
   params[[2]] <- list(X = X, Y = Y, ncomp = 1)
   params[[3]] <- list(X = X, Y = Y, ncomp = A)
   params[[4]] <- list(X = X, Y = Y, ncomp = A, cv = 10)
   params[[5]] <- list(X = X, Y = Y, ncomp = A, cv = 10, scale = TRUE)
   params[[6]] <- list(X = X, Y = Y, ncomp = A, cv = list("ven", 4))
   params[[7]] <- list(X = X, Y = Y, ncomp = A, cv = list("ven", 4), scale = TRUE)
   params[[8]] <- list(X = X, Y = Y, ncomp = A, cv = list("loo"), scale = TRUE)
   params[[9]] <- list(X = X, Y = Y, ncomp = A, cv = list("loo"))

   for (i in seq_along(params)) {
      expect_silent(Xpv <- do.call(pcvpcr, params[[i]]))
      expect_equal(nrow(Xpv), nrow(X))
      expect_equal(ncol(Xpv), ncol(X))

      D <- attr(Xpv, "D")
      expect_false(is.null(D))
      expect_equal(ncol(D), if (is.null(params[[i]]$ncomp)) 30 else params[[i]]$ncomp)
   }
})

test_that("- pcvpcr() works well for Corn data.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv = rep(seq_len(4), length.out = nrow(X))[order(Y)]

   Xpv <- pcvpcr(X, Y, 20, center = TRUE, scale = FALSE, cv = cv)
   D <- attr(Xpv, "D")

   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))
   expect_false(is.null(D))
   expect_equal(ncol(D), 20)
   expect_equal(nrow(D), 4)

   mX = apply(X, 2, mean)
   Xmc = scale(X, center = mX, scale = FALSE)
   Xpvmc = scale(Xpv, center = mX, scale = FALSE)
   Ymc = scale(Y, center = TRUE, scale = FALSE)

   P = svd(Xmc)$v[, 1:20]
   T = Xmc %*% P;

   Tpv = Xpvmc %*% P;

   rmse <- rep(0, 20)
   for (a in 1:20) {
      C <- solve(crossprod(T[, 1:a, drop = FALSE])) %*% crossprod(T[, 1:a, drop = FALSE], Ymc)
      Ypv <- Tpv[, 1:a, drop = FALSE] %*% C;
      Ey <- Ymc - Ypv
      rmse[a] <- sqrt(sum(Ey^2) / nrow(X))
   }

   expect_equivalent(round(rmse, 3), c(0.363, 0.331, 0.298, 0.298, 0.271, 0.251, 0.213, 0.213, 0.214, 0.216, 0.216, 0.209, 0.190, 0.181, 0.183, 0.187, 0.187, 0.183, 0.184, 0.193))
})

context("Tests for 'pcvpls()':")

test_that("- pcvpls() works well for random data.", {

   I <- 100
   J <- 50
   A <- 30

   set.seed(42)
   X <- matrix(rnorm(I * J), I, J)
   Y <- X %*% runif(J)

   params <- list()
   params[[1]] <- list(X = X, Y = Y)
   params[[2]] <- list(X = X, Y = Y, ncomp = 1)
   params[[3]] <- list(X = X, Y = Y, ncomp = A)
   params[[4]] <- list(X = X, Y = Y, ncomp = A, cv = 10)
   params[[5]] <- list(X = X, Y = Y, ncomp = A, cv = 10, scale = TRUE)
   params[[6]] <- list(X = X, Y = Y, ncomp = A, cv = list("ven", 4))
   params[[7]] <- list(X = X, Y = Y, ncomp = A, cv = list("ven", 4), scale = TRUE)
   params[[8]] <- list(X = X, Y = Y, ncomp = A, cv = list("loo"), scale = TRUE)
   params[[9]] <- list(X = X, Y = Y, ncomp = A, cv = list("loo"))

   for (i in seq_along(params)) {
      expect_silent(Xpv <- do.call(pcvpls, params[[i]]))
      expect_equal(nrow(Xpv), nrow(X))
      expect_equal(ncol(Xpv), ncol(X))

      D <- attr(Xpv, "D")
      expect_false(is.null(D))
      expect_equal(ncol(D), if (is.null(params[[i]]$ncomp)) 30 else params[[i]]$ncomp)
   }
})

test_that("- pcvpls() works well for Corn data with global CV scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv = rep(seq_len(4), length.out = nrow(X))[order(Y)]
   Xpv <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv)
   D <- attr(Xpv, "D")

   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))
   expect_false(is.null(D))
   expect_equal(ncol(D), 20)
   expect_equal(nrow(D), 4)
   expect_true(min(D) > 0.0)
   expect_true(max(D) < 1.5)
   expect_true(all(D[, 4] > 1.1))
   expect_silent(plotD(Xpv))
})

test_that("- pcvpls() works well for Corn data with local CV scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv = rep(seq_len(4), length.out = nrow(X))[order(Y)]
   Xpv <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "local")
   D <- attr(Xpv, "D")

   mdatools::mdaplot(Xpv, type = "l")
   plotD(Xpv)
   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))
   expect_false(is.null(D))
   expect_equal(ncol(D), 20)
   expect_equal(nrow(D), 4)
   expect_true(min(D) > 0.0)
   expect_true(max(D) < 1.5)
   expect_true(all(D[, 4] > 1.1))
   expect_silent(plotD(Xpv))
})

test_that("- pcvpls() compare local and global CV scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   cv = list("ven", 10)
   Xpv1 <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "global")
   Xpv2 <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "local")

   par(mfrow = c(2, 2))
   mdatools::mdaplot(Xpv1, type = "l", cgroup = Y[, 1])
   mdatools::mdaplot(Xpv2, type = "l", cgroup = Y[, 1])
   plotD(Xpv1)
   plotD(Xpv2)

   m11 = mdatools::pls(X, Y, 20, cv = cv, cv.scope = "global")
   m21 = mdatools::pls(X, Y, 20, cv = cv, cv.scope = "local")
   m12 = mdatools::pls(X, Y, 20, x.test = Xpv1, y.test = Y, cv.scope = "global")
   m22 = mdatools::pls(X, Y, 20, x.test = Xpv2, y.test = Y, cv.scope = "local")

   expect_equivalent(m11$res$cv$rmse, m12$res$test$rmse)
   expect_equivalent(m21$res$cv$rmse, m22$res$test$rmse, tolerance = 0.01)

})
