####################################
# Tests for PCV for PCA/SIMCA      #
####################################

# global settings

## directory to keep CSV files with reference values
caseDir <- "../../../.tests/pcvpca/"

setup({
   pdf(file = tempfile("pcv-test-pca-", fileext = ".pdf"))
})

teardown({
   dev.off()
})

#' Show distance plot for the results
plotdistance <- function(rcl, rcv, rpv, a, main = "Distance plot") {

   q <- rcl$Q[, a]
   qpv <- rpv$Q[, a]
   qcv <- rcv$Q[, a]
   q0 <- mean(q)

   h <- rcl$H[, a]
   hpv <- rpv$H[, a]
   hcv <- rcv$H[, a]
   h0 <- mean(h)

   plot(
      c(h/h0, hcv/h0, hpv/h0), c(q/q0, qcv/q0, qpv/q0),
      col = rep(c("blue", "red", "black"), each = length(h)),
      pch = rep(c(5, 1, 4), each = length(h)),
      xlab = "Score distance, h/h0",
      ylab = "Orthogonal distance, q/q0",
      main = main
   )

}

#' Create global PCA model and apply it to Xpv set
pcapv <- function(X, Xpv, ncomp, center = TRUE, scale = FALSE) {

   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale)  apply(X, 2, sd) else rep(1, ncol(X))

   Xs <- scale(X, center = mX, scale = sX)
   Xpvs <- scale(Xpv, center = mX, scale = sX)

   # global model
   m <- svd(Xs, nv = ncomp, nu = ncomp)
   s <- m$d[seq_len(ncomp)] / sqrt(nrow(X) - 1)
   P <- m$v

   # predictions for PV-set
   Tpv <- Xpvs %*% P;
   Upv <- Tpv %*% diag(1 / s, ncomp, ncomp);

   Qpv <- matrix(0, nrow(X), ncomp)
   Hpv <- matrix(0, nrow(X), ncomp)

   for (a in seq_len(ncomp)) {

      aind <- seq_len(a)
      P.a <- P[, aind, drop = FALSE]
      Tpv.a <- Tpv[, aind, drop = FALSE]
      Upv.a <- Upv[, aind, drop = FALSE]

      Epv.a <- Xpvs - tcrossprod(Tpv.a, P.a)
      Qpv[, a] <- rowSums(Epv.a^2)
      Hpv[, a] <- rowSums(Upv.a^2)
   }

   return (list(Q = Qpv, H = Hpv))
}

#' Do cross-validation of PCA model with global scope
pcacvglobal <- function(X, ncomp, cv, center = TRUE, scale = FALSE) {

   cvind <- pcvcrossval(cv, nrow(X), seq_len(nrow(X)))

   # scale data globally
   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale)  apply(X, 2, sd) else rep(1, ncol(X))
   Xs <- scale(X, center = mX, scale = sX)

   # fit global model
   m <- svd(Xs, nv = ncomp, nu = ncomp)
   s <- m$d[seq_len(ncomp)] / sqrt(nrow(X) - 1)

   nseg <- max(cvind)
   Qcv <- matrix(0, nrow(X), ncomp)
   Hcv <- matrix(0, nrow(X), ncomp)

   for (k in seq_len(nseg)) {

      # get indices for local sets
      indc <- cvind != k
      indk <- cvind == k

      # split data to local calibration and validation sets
      Xc <- Xs[indc, , drop = FALSE]
      Xk <- Xs[indk, , drop = FALSE]

      # fit the local model
      mk <- svd(Xc, nv = ncomp, nu = ncomp)
      Pk <- mk$v

      # make predictions for local validation set
      Tk <- Xk %*% Pk;
      Uk <- Tk %*% diag(1 / s, ncomp, ncomp)

      # compute distances
      for (a in seq_len(ncomp)) {

         aind <- seq_len(a)
         Tk.a <- Tk[, aind, drop = FALSE]
         Uk.a <- Uk[, aind, drop = FALSE]
         Pk.a <- Pk[, aind, drop = FALSE]

         Ek.a <- Xk - tcrossprod(Tk.a, Pk.a)
         Qcv[indk, a] <- rowSums(Ek.a^2)
         Hcv[indk, a] <- rowSums(Uk.a^2)
      }
   }

   return (list(Q = Qcv, H = Hcv))
}

#' Do cross-validation of PCA model with local scope
pcacvlocal <- function(X, ncomp, cv, center = TRUE, scale = FALSE) {

   cvind <- pcvcrossval(cv, nrow(X), seq_len(nrow(X)))

   # scale data globally
   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale)  apply(X, 2, sd) else rep(1, ncol(X))
   Xs <- scale(X, center = mX, scale = sX)

   # fit global model to gen singular values
   m <- svd(Xs, nv = ncomp, nu = ncomp)
   s <- m$d[seq_len(ncomp)] / sqrt(nrow(X) - 1)

   nseg <- max(cvind)
   Qcv <- matrix(0, nrow(X), ncomp)
   Hcv <- matrix(0, nrow(X), ncomp)

   for (k in seq_len(nseg)) {

      # get indices for local sets
      indc <- cvind != k
      indk <- cvind == k

      # split data to local calibration and validation sets
      Xc <- X[indc, , drop = FALSE]
      Xk <- X[indk, , drop = FALSE]

      # local [re]-scaling
      mX <- if (center) apply(Xc, 2, mean) else rep(0, ncol(Xc))
      sX <- if (scale)  apply(Xc, 2, sd) else rep(1, ncol(Xc))
      Xc <- scale(Xc, center = mX, scale = sX)
      Xk <- scale(Xk, center = mX, scale = sX)

      # fit local model
      mk <- svd(Xc, nv = ncomp, nu = ncomp)
      Pk <- mk$v

      # make predictions for local validation set
      Tk <- Xk %*% Pk;
      Uk <- Tk %*% diag(1 / s, ncomp, ncomp)

      # compute distances
      for (a in seq_len(ncomp)) {

         aind <- seq_len(a)
         Tk.a <- Tk[, aind, drop = FALSE]
         Uk.a <- Uk[, aind, drop = FALSE]
         Pk.a <- Pk[, aind, drop = FALSE]

         Ek.a <- Xk - tcrossprod(Tk.a, Pk.a)
         Qcv[indk, a] <- rowSums(Ek.a^2)
         Hcv[indk, a] <- rowSums(Uk.a^2)
      }
   }

   return (list(Q = Qcv, H = Hcv))
}

#' Save results as reference values
savereferences <- function(rpv.g, rpv.l, scale, ncomp, cv) {
   cvText <- if (length(cv) == 2) paste0(cv[[1]], cv[[2]]) else cv[[1]]
   caseSuffix <- sprintf("-%d-%s-%s.csv", ncomp, scale, cvText)

   write.table(rpv.g$Q, file = paste0(caseDir, "Qpvg", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(rpv.g$H, file = paste0(caseDir, "Hpvg", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(rpv.l$Q, file = paste0(caseDir, "Qpvl", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(rpv.l$H, file = paste0(caseDir, "Hpvl", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
}

#' Run tests for different combinations of parameters
runtests <- function(X, save.res = FALSE) {

   # test cases
   cv_cases <- list(list("ven", 4), list("ven", 10), list("loo"), list("rand", 10))
   ncomp_cases <- if (ncol(X) < 40) c(1, 4, 8, 12) else c(1, 10, 20, 30)
   scale_cases <- c(TRUE, FALSE)
   center <- TRUE

   # run tests
   cases <- expand.grid(cv = cv_cases, ncomp = ncomp_cases, scale = scale_cases)
   for (i in seq_len(nrow(cases))) {

      # get settings for current case
      cv <- cases[i, "cv"][[1]]
      ncomp <- cases[i, "ncomp"]
      scale <- cases[i, "scale"]

      # tets results for global scope
      set.seed(42)
      Xpv.g <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "global")

      set.seed(42)
      rcv.g <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv.g <- pcapv(X, Xpv.g, ncomp = ncomp, center = center, scale = scale)

      expect_equal(rcv.g$H, rpv.g$H, check.attributes = FALSE)
      expect_equal(rcv.g$Q, rpv.g$Q, check.attributes = FALSE)

      # test results for local scope
      set.seed(42)
      Xpv.l <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")

      set.seed(42)
      rcv.l <- pcacvlocal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv.l <- pcapv(X, Xpv.l, ncomp = ncomp, center = center, scale = scale)

      expect_equal(rcv.l$H, rpv.l$H, check.attributes = FALSE)
      expect_equal(rcv.l$Q, rpv.l$Q, check.attributes = FALSE)

      # save outcomes as reference
      if (save.res && cv[[1]] != "rand") {
         savereferences(rpv.g, rpv.l, scale, ncomp, cv)
      }
   }
}


test_that("automatic tests for People data are passed.", {

   load("People.RData")
   runtests(people)

})

test_that("automatic tests for Corn data are passed.", {

   data(corn)
   X <- corn$spectra

   # run automatic tests and save results as reference values
   runtests(X, save.res = TRUE)
})

test_that("manual tests for Corn data are passed.", {

   data(corn)
   X <- corn$spectra

   # manual tests
   ncomp <- 30
   center <- TRUE
   scale <- FALSE
   cv <- list("ven", 4)

   Xpv <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv)

   rcl <- pcapv(X, X, ncomp = ncomp, center = center, scale = scale)
   rpv <- pcapv(X, Xpv, ncomp = ncomp, center = center, scale = scale)
   rcv <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)

   a <- 2;

   q <- rcl$Q[, a]
   qpv <- rpv$Q[, a]
   h <- rcl$H[, a]
   hpv <- rpv$H[, a]

   expect_equal(sum(q/mean(q) > 4), 4)
   expect_equal(sum(qpv/mean(q) > 4), 4)
   expect_equal(sum(h/mean(h) > 5), 2)
   expect_equal(sum(h/mean(h) > 5), 2)

   a <- 20;

   q <- rcl$Q[, a]
   qpv <- rpv$Q[, a]
   h <- rcl$H[, a]
   hpv <- rpv$H[, a]

   expect_equal(sum(q/mean(q) > 4), 0)
   expect_equal(sum(qpv/mean(q) > 4), 5)
   expect_equal(sum(h/mean(h) > 2), 3)
   expect_equal(sum(hpv/mean(h) > 2), 9)

})

test_that("visual tests for Corn data looks good.", {

   data(corn)
   X <- corn$spectra

   ncomp <- 30
   center <- TRUE
   scale <- TRUE
   cv <- list("ven", 20)

   # global scope
   Xpv.g <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv)
   rcv.g <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rpv.g <- pcapv(X, Xpv.g, ncomp = ncomp, center = center, scale = scale)
   rcl.g <- pcapv(X, X, ncomp = ncomp, center = center, scale = scale)

   # local scope
   Xpv.l <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")
   rcv.l <- pcacvlocal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rpv.l <- pcapv(X, Xpv.l, ncomp = ncomp, center = center, scale = scale)
   rcl.l <- pcapv(X, X, ncomp = ncomp, center = center, scale = scale)

   # distance plots
   par(mfrow = c(2, 2))
   plotdistance(rcl.g, rcv.g, rpv.g,  4, main = "Global, a = 4")
   plotdistance(rcl.g, rcv.g, rpv.g, 20, main = "Global, a = 20")
   plotdistance(rcl.l, rcv.l, rpv.l,  4, main = "Local, a = 4")
   plotdistance(rcl.l, rcv.l, rpv.l, 20, main = "Local, a = 20")

})

test_that("edge cases: constant column, small data, invalid cv.scope.", {

   data(corn)
   X <- corn$spectra

   # data with a constant column appended — tests epsilon protection on std dev
   Xc <- cbind(X, rep(42, nrow(X)))
   nc <- ncol(Xc)
   set.seed(42)
   Xpv <- pcvpca(Xc, ncomp = 5, scale = TRUE, cv = list("ven", 4))
   expect_equal(nrow(Xpv), nrow(Xc))
   expect_equal(ncol(Xpv), ncol(Xc))
   expect_false(any(is.nan(Xpv)))

   # the constant column must be reproduced as is, not dropped or zeroed: after centering it
   # is exactly zero, contributes nothing to the model, and uncentering restores the constant
   expect_equal(Xpv[, nc], rep(42, nrow(Xc)), ignore_attr = TRUE)

   # same must hold for the local scope
   set.seed(42)
   Xpv <- pcvpca(Xc, ncomp = 5, scale = TRUE, cv = list("ven", 4), cv.scope = "local")
   expect_equal(Xpv[, nc], rep(42, nrow(Xc)), ignore_attr = TRUE)

   # without centering the constant column carries variance for the uncentered decomposition,
   # so it is a part of the model — here we only require finite values of a correct size
   set.seed(42)
   Xpv <- pcvpca(Xc, ncomp = 5, center = FALSE, scale = TRUE, cv = list("ven", 4))
   expect_equal(dim(Xpv), dim(Xc))
   expect_true(all(is.finite(Xpv)))

   # small dataset (6 rows) — tests ncomp auto-reduction
   Xs <- X[1:6, ]
   set.seed(42)
   Xpv <- pcvpca(Xs, ncomp = 20, cv = list("ven", 3))
   expect_equal(nrow(Xpv), 6)
   expect_false(any(is.nan(Xpv)))

   # invalid cv.scope should error
   expect_error(pcvpca(X, ncomp = 5, cv.scope = "wrong"), "should be one of")
})

test_that("uneven split: number of components is limited by the largest segment.", {
   set.seed(1)
   X <- matrix(rnorm(20 * 10), 20, 10)
   cv <- c(rep(1, 16), rep(2, 4))

   # the largest segment leaves 4 rows, so at most 3 components are possible
   set.seed(2)
   Xpv9 <- pcvpca(X, 9, cv = cv)
   set.seed(2)
   Xpv3 <- pcvpca(X, 3, cv = cv)
   expect_true(all(is.finite(Xpv9)))
   expect_equal(Xpv9, Xpv3)

   # orthogonal distances of the PV-set on the global model equal cross-validation
   Xc <- scale(X, center = TRUE, scale = FALSE)
   P <- svd(Xc, nu = 0, nv = 3)$v
   Xpvc <- sweep(Xpv3, 2, attr(Xc, "scaled:center"))
   qpv <- rowSums((Xpvc - Xpvc %*% P %*% t(P))^2)
   qcv <- numeric(20)
   for (k in 1:2) {
      Pk <- svd(Xc[cv != k, ], nu = 0, nv = 3)$v
      Xk <- Xc[cv == k, , drop = FALSE]
      qcv[cv == k] <- rowSums((Xk - Xk %*% Pk %*% t(Pk))^2)
   }
   expect_equal(qpv, qcv)

   # several repetitions: only the first one is used
   expect_warning(pcvpca(X, 3, cv = list("rand", 4, 3)), "first one")

   # local calibration set with one row: no component possible
   expect_error(pcvpca(X, 1, cv = c(rep(1, 19), 2)), "at least 2 rows")
})
