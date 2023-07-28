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

      expect_equivalent(rcv.g$H, rpv.g$H)
      expect_equivalent(rcv.g$Q, rpv.g$Q)

      # test results for local scope
      set.seed(42)
      Xpv.l <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")

      set.seed(42)
      rcv.l <- pcacvlocal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv.l <- pcapv(X, Xpv.l, ncomp = ncomp, center = center, scale = scale)

      expect_equivalent(rcv.l$H, rpv.l$H)
      expect_equivalent(rcv.l$Q, rpv.l$Q)

      # save outcomes as reference
      if (save.res && cv[[1]] != "rand") {
         savereferences(rpv.g, rpv.l, scale, ncomp, cv)
      }
   }
}


context("Tests for 'pcvpca()':")

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
