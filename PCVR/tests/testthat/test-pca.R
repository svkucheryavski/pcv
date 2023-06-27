####################################
# Tests for PCV for PCA/SIMCA      #
####################################

setup({
   pdf(file = tempfile("pcv-test-pca-", fileext = ".pdf"))
})

teardown({
   dev.off()
})

#' Create global PCA model and apply it to Xpv set
pcapv <- function(X, Xpv, ncomp = min(nrow(X) - 1, ncol(X), 30), center = TRUE, scale = FALSE) {

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
pcacvglobal <- function(X, ncomp = min(nrow(X) - 1, ncol(X), 30), cv = list("ven", 4), center = TRUE, scale = FALSE) {

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
pcacvlocal <- function(X, ncomp = min(nrow(X) - 1, ncol(X), 30), cv = list("ven", 4), center = TRUE, scale = FALSE) {

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
      Xc <- Xs[indc, , drop = FALSE]
      Xk <- Xs[indk, , drop = FALSE]

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


context("Tests for 'pcvpca()':")

test_that("- pcvpca() works well for random data - simple test.", {

   I <- 100
   J <- 50
   X <- matrix(rnorm(I * J), I, J)

   ncomp <- 30
   cv <- list("ven", 4)
   center <- TRUE
   scale <- FALSE

   Xpv <- pcvpca(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rcv <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rpv <- pcapv(X, Xpv, ncomp = ncomp, center = center, scale = scale)

   expect_equivalent(rcv$H, rpv$H)
   expect_equivalent(rcv$Q, rpv$Q)
})

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

      X <- params[[i]]$X
      ncomp <- if (is.null(params[[i]]$ncomp)) min(nrow(X) - 1, ncol(X), 30) else params[[i]]$ncomp
      cv <- if (is.null(params[[i]]$cv)) list("ven", 4) else params[[i]]$cv
      center <- if (is.null(params[[i]]$center)) TRUE else params[[i]]$center
      scale <- if (is.null(params[[i]]$scale)) FALSE else params[[i]]$scale

      # this is needed for reproducibility if cv is random
      set.seed(42)
      expect_silent(Xpv <- do.call(pcvpca, params[[i]]))
      expect_equal(nrow(Xpv), nrow(X))
      expect_equal(ncol(Xpv), ncol(X))
      expect_true(ks.test(X, Xpv)$p.value > 0.01)


      # this is needed for reproducibility if cv is random
      set.seed(42)
      rcv <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv <- pcapv(X, Xpv, ncomp = ncomp, center = center, scale = scale)

      expect_equivalent(rcv$H, rpv$H)
      expect_equivalent(rcv$Q, rpv$Q)
   }

})

test_that("- pcvpca() works well for Corn data.", {
   data(corn)
   X <- corn$spectra

   Xpv <- pcvpca(X, 30, center = TRUE, scale = FALSE, cv = list("ven", 4))
   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))

   mX <- apply(X, 2, mean)
   Xmc <- scale(X, center = mX, scale = FALSE)
   Xpvmc <- scale(Xpv, center = mX, scale = FALSE)

   P <- svd(Xmc)$v[, 1:20]
   T <- Xmc %*% P;
   lambda <- colSums(T^2 / (nrow(X) - 1))

   Tpv <- Xpvmc %*% P;

   # check the distances (reproducing example from the paper)
   a <- 2;
   E <- Xmc - tcrossprod(T[, 1:a], P[, 1:a]);
   Epv <- Xpvmc - tcrossprod(Tpv[, 1:a], P[, 1:a]);
   q <- rowSums(E^2)
   qpv <- rowSums(Epv^2)
   h <- rowSums(T[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   hpv <- rowSums(Tpv[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   expect_equal(sum(q/mean(q) > 4), 4)
   expect_equal(sum(qpv/mean(q) > 4), 4)
   expect_equal(sum(h/mean(h) > 5), 2)
   expect_equal(sum(h/mean(h) > 5), 2)

   a <- 20;
   E <- Xmc - tcrossprod(T[, 1:a], P[, 1:a]);
   Epv <- Xpvmc - tcrossprod(Tpv[, 1:a], P[, 1:a]);
   q <- rowSums(E^2)
   qpv <- rowSums(Epv^2)
   h <- rowSums(T[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))
   hpv <- rowSums(Tpv[, 1:a]^2 %*% diag(1/lambda[1:a], a, a))

   expect_equal(sum(q/mean(q) > 4), 0)
   expect_equal(sum(qpv/mean(q) > 4), 5)
   expect_equal(sum(h/mean(h) > 2), 3)
   expect_equal(sum(hpv/mean(h) > 2), 9)

   # compare with manual computation
   set.seed(42)
   rcv <- pcacvglobal(X, ncomp = 30, cv = list("ven", 4), center = TRUE, scale = FALSE)
   rpv <- pcapv(X, Xpv, ncomp = 30, center = TRUE, scale = FALSE)

   expect_equivalent(rcv$H, rpv$H)
   expect_equivalent(rcv$Q, rpv$Q)

})

test_that("- pcvpca() works well for Corn data and both global and local scope.", {
   data(corn)
   X <- corn$spectra

   ncomp <- 30
   cv <- list("ven", 10)
   center <- TRUE
   scale <- FALSE

   # generate PV-sets
   Xpv.l <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")
   Xpv.g <- pcvpca(X, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "global")

   # get CCV results
   rcv.l <- pcacvlocal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rcv.g <- pcacvglobal(X, ncomp = ncomp, cv = cv, center = center, scale = scale)

   # get PCV results
   rpv.l <- pcapv(X, Xpv.l, ncomp = ncomp, center = center, scale = scale)
   rpv.g <- pcapv(X, Xpv.g, ncomp = ncomp, center = center, scale = scale)

   # here we expect that the difference between PV and CCV results will not be
   # more than 5% of maximum distance for the same component
   expect_lte(max(sweep(abs(rcv.l$H - rpv.l$H), 2, apply(rcv.l$H, 2, max), "/")),  0.05)
   expect_lte(max(sweep(abs(rcv.l$Q - rpv.l$Q), 2, apply(rcv.l$Q, 2, max), "/")),  0.05)

   # compare global PV-results with CCV
   expect_equivalent(rcv.g$H, rpv.g$H)
   expect_equivalent(rcv.g$Q, rpv.g$Q)

   # do PCA and get PV results from there
   m <- mdatools::pca(X, 20, center = TRUE, scale = FALSE)
   r.l <- predict(m, Xpv.l)
   r.g <- predict(m, Xpv.g)

   # show plots
   par(mfrow = c(3, 2))

   # distance plots for local results
   mdatools::plotResiduals(m, pch = 1, log = FALSE, norm = FALSE, res = list(cal = m$calres, pv = r.l), ncomp =  2)
   points(rcv.l$H[, 2], rcv.l$Q[, 2], pch = 4)
   mdatools::plotResiduals(m, pch = 1, log = FALSE, norm = FALSE, res = list(cal = m$calres, pv = r.l), ncomp = 10)
   points(rcv.l$H[, 10], rcv.l$Q[, 10], pch = 4)

   # distance plots for global results
   mdatools::plotResiduals(m, pch = 1, log = FALSE, norm = FALSE, res = list(cal = m$calres, pv = r.g), ncomp =  2)
   points(rcv.g$H[, 2], rcv.g$Q[, 2], pch = 4)
   mdatools::plotResiduals(m, pch = 1, log = FALSE, norm = FALSE, res = list(cal = m$calres, pv = r.g), ncomp = 10)
   points(rcv.g$H[, 10], rcv.g$Q[, 10], pch = 4)

   # line plots with PV-sets
   mdatools::mdaplot(mdatools::prep.autoscale(Xpv.l), type = "l")
   mdatools::mdaplot(mdatools::prep.autoscale(Xpv.g), type = "l")


})

