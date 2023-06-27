####################################
# Tests for PCV for PLS            #
####################################

setup({
   pdf(file = tempfile("pcv-test-pls-", fileext = ".pdf"))
})

teardown({
   dev.off()
})


#' Create global PCR model and apply it to Xpv set
plspv <- function(X, Y, Xpv, ncomp, center = TRUE, scale = FALSE) {

   if (is.null(dim(Y))) {
      dim(Y) <- c(nrow(X), 1)
   }

   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale)  apply(X, 2, sd) else rep(1, ncol(X))
   mY <- if (center) apply(Y, 2, mean) else 0
   sY <- if (center) apply(Y, 2, sd) else 1

   Xs <- scale(X, center = mX, scale = sX)
   Ys <- scale(Y, center = mY, scale = sY)
   Xpvs <- scale(Xpv, center = mX, scale = sX)

   # global model
   m <- simpls(Xs, Ys, ncomp)
   T <- Xs %*% m$R;

   # predictions for PV-set
   Tpv <- Xpvs %*% m$R;
   Ypv <- matrix(0, nrow(Y), ncomp)

   for (a in 1:ncomp) {
      Ypv[, a] <- tcrossprod(Tpv[, 1:a, drop = FALSE], m$C[, 1:a, drop = FALSE])
   }

   Ypv <- Ypv * sY + mY
   E <- apply(Ypv, 2, function(v) v - Y)

   return (list(Yp = Ypv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

#' Do cross-validation of PCR model with global scope
plscvglobal <- function(X, Y, ncomp, cv, center = TRUE, scale = FALSE) {

   if (is.null(dim(Y))) {
      dim(Y) <- c(nrow(X), 1)
   }

   cvind <- pcvcrossval(cv, nrow(X), Y)

   mX <- if (center) apply(X, 2, mean) else rep(0, ncol(X))
   sX <- if (scale)  apply(X, 2, sd) else rep(1, ncol(X))
   mY <- if (center) apply(Y, 2, mean) else 0
   sY <- if (center) apply(Y, 2, sd) else 1

   Xs <- scale(X, center = mX, scale = sX)
   Ys <- scale(Y, center = mY, scale = sY)

   nSeg <- max(cvind)
   Ycv <- matrix(0, nrow(Y), ncomp)

   for (i in unique(cvind)) {

      indc <- cvind != i
      indk <- cvind == i

      Xc <- Xs[indc, , drop = FALSE]
      Yc <- Ys[indc, , drop = FALSE]
      Xk <- Xs[indk, , drop = FALSE]

      mk <- simpls(Xc, Yc, ncomp)
      Tc <- Xc %*% mk$R;
      Tk <- Xk %*% mk$R;

      for (a in 1:ncomp) {
         Ycv[indk, a] <- tcrossprod(Tk[, 1:a, drop = FALSE], mk$C[, 1:a, drop = FALSE])
      }
   }

   Ycv <- Ycv * sY + mY
   E <- apply(Ycv, 2, function(v) v - Y)

   return (list(Yp = Ycv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

#' Do cross-validation of PCR model with local scope
plscvlocal <- function(X, Y, ncomp, cv, center = TRUE, scale = FALSE) {

   if (is.null(dim(Y))) {
      dim(Y) <- c(nrow(X), 1)
   }

   cvind <- pcvcrossval(cv, nrow(X), Y)

   nSeg <- max(cvind)
   Ycv <- matrix(0, nrow(Y), ncomp)
   for (i in unique(cvind)) {

      indc <- cvind != i
      indk <- cvind == i

      Xc <- X[indc, , drop = FALSE]
      Yc <- Y[indc, , drop = FALSE]
      Xk <- X[indk, , drop = FALSE]

      mX <- if (center) apply(Xc, 2, mean) else rep(0, ncol(Xc))
      sX <- if (scale)  apply(Xc, 2, sd) else rep(1, ncol(Xc))
      mY <- if (center) apply(Yc, 2, mean) else 0
      sY <- if (center) apply(Yc, 2, sd) else 1

      Xcs <- scale(Xc, center = mX, scale = sX)
      Ycs <- scale(Yc, center = mY, scale = sY)
      Xks <- scale(Xk, center = mX, scale = sX)


      mk <- simpls(Xcs, Ycs, ncomp)
      Tc <- Xcs %*% mk$R;
      Tk <- Xks %*% mk$R;

      for (a in 1:ncomp) {
         Ycv[indk, a] <- tcrossprod(Tk[, 1:a, drop = FALSE], mk$C[, 1:a, drop = FALSE]);
      }

      Ycv[indk, ] <- Ycv[indk, ] * sY + mY;
   }

   E <- apply(Ycv, 2, function(v) v - Y)

   return (list(Yp = Ycv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}


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

      Rk <- attr(Xpv, "Rk")
      expect_false(is.null(Rk))

      R <- attr(Xpv, "R")
      expect_false(is.null(R))

   }
})

test_that("- pcvpls() works well for Corn data with global CV scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv <- rep(seq_len(4), length.out = nrow(X))[order(Y)]
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

   cv <- list("ven", 10)
   Xpv.g <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "global")
   Xpv.l <- pcvpls(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "local")

   par(mfrow = c(2, 2))
   matplot(t(scale(Xpv.g, center = TRUE, scale = FALSE)), type = "l", lty = 1, main = "Global")
   matplot(t(scale(Xpv.l, center = TRUE, scale = FALSE)), type = "l", lty = 1, main = "Local")
   plotD(Xpv.g, main = "Global")
   plotD(Xpv.l, main = "Local")

   rcv.g <- plscvglobal(X, Y, 20, cv = cv)
   rcv.l <- plscvlocal(X, Y, 20, cv = cv)
   rpv.g <- plspv(X, Y, Xpv.g, 20)
   rpv.l <- plspv(X, Y, Xpv.l, 20)

   expect_equivalent(rcv.g$RMSE, rpv.g$RMSE)
   expect_equivalent(rcv.l$RMSE, rpv.l$RMSE, tolerance = 0.01)

})
