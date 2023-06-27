####################################
# Tests for PCV for PCR            #
####################################

setup({
   pdf(file = tempfile("pcv-test-pcr-", fileext = ".pdf"))
})

teardown({
   dev.off()
})

#' Create global PCR model and apply it to Xpv set
pcrpv <- function(X, Y, Xpv, ncomp, center = TRUE, scale = FALSE) {

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
   P <- svd(Xs)$v[, 1:ncomp]
   T <- Xs %*% P;

   # predictions for PV-set
   Tpv <- Xpvs %*% P;
   Ypv <- matrix(0, nrow(Y), ncomp)

   for (a in 1:ncomp) {
      C <- solve(crossprod(T[, 1:a, drop = FALSE])) %*% crossprod(T[, 1:a, drop = FALSE], Ys)
      Ypv[, a] <- Tpv[, 1:a, drop = FALSE] %*% C;
   }

   Ypv <- Ypv * sY + mY
   E <- apply(Ypv, 2, function(v) v - Y)

   return (list(Yp = Ypv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

#' Do cross-validation of PCR model with global scope
pcrcvglobal <- function(X, Y, ncomp, cv, center = TRUE, scale = FALSE) {

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

      Pk <- svd(Xc)$v[, 1:ncomp]
      Tc <- Xc %*% Pk;
      Tk <- Xk %*% Pk;

      for (a in 1:ncomp) {
         C <- solve(crossprod(Tc[, 1:a, drop = FALSE])) %*% crossprod(Tc[, 1:a, drop = FALSE], Yc)
         Ycv[indk, a] <- Tk[, 1:a, drop = FALSE] %*% C;
      }
   }

   Ycv <- Ycv * sY + mY
   E <- apply(Ycv, 2, function(v) v - Y)

   return (list(Yp = Ycv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

#' Do cross-validation of PCR model with local scope
pcrcvlocal <- function(X, Y, ncomp, cv, center = TRUE, scale = FALSE) {

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


      Pk = svd(Xcs)$v[, 1:ncomp]
      Tc = Xcs %*% Pk;
      Tk = Xks %*% Pk;

      for (a in 1:ncomp) {
         C <- solve(crossprod(Tc[, 1:a, drop = FALSE])) %*% crossprod(Tc[, 1:a, drop = FALSE], Ycs)
         Ycv[indk, a] <- Tk[, 1:a, drop = FALSE] %*% C;
      }

      Ycv[indk, ] <- Ycv[indk, ] * sY + mY;
   }

   E <- apply(Ycv, 2, function(v) v - Y)

   return (list(Yp = Ycv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

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

      X <- params[[i]]$X
      Y <- params[[i]]$Y
      ncomp <- if (is.null(params[[i]]$ncomp)) min(nrow(X) - 1, ncol(X), 30) else params[[i]]$ncomp
      cv <- if (is.null(params[[i]]$cv)) list("ven", 4) else params[[i]]$cv
      center <- if (is.null(params[[i]]$center)) TRUE else params[[i]]$center
      scale <- if (is.null(params[[i]]$scale)) FALSE else params[[i]]$scale

      # this is needed for reproducibility if cv is random
      set.seed(42)
      expect_silent(Xpv <- do.call(pcvpcr, params[[i]]))
      expect_equal(nrow(Xpv), nrow(X))
      expect_equal(ncol(Xpv), ncol(X))

      D <- attr(Xpv, "D")
      expect_false(is.null(D))
      expect_equal(ncol(D), ncomp)

      # this is needed for reproducibility if cv is random
      set.seed(42)
      rcv <- pcrcvglobal(X, Y, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv <- pcrpv(X, Y, Xpv, ncomp = ncomp, center = center, scale = scale)

      expect_equivalent(rcv$Yp, rpv$Yp)
      expect_equivalent(rcv$RMSE, rpv$RMSE)
   }
})

test_that("- pcvpcr() works well for Corn data.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv <- rep(seq_len(4), length.out = nrow(X))[order(Y)]

   Xpv <- pcvpcr(X, Y, 20, center = TRUE, scale = FALSE, cv = cv)
   D <- attr(Xpv, "D")

   expect_equal(nrow(Xpv), nrow(X))
   expect_equal(ncol(Xpv), ncol(X))
   expect_false(is.null(D))
   expect_equal(ncol(D), 20)
   expect_equal(nrow(D), 4)

   rcv <- pcrcvglobal(X, Y, ncomp = 20, cv = cv, center = TRUE, scale = FALSE)
   rpv <- pcrpv(X, Y, Xpv, ncomp = 20, center = TRUE, scale = FALSE)

   expect_equivalent(rcv$Yp, rpv$Yp)
   expect_equivalent(rcv$RMSE, rpv$RMSE)
   expect_equivalent(round(rpv$RMSE, 3),
      c(0.363, 0.331, 0.298, 0.298, 0.271, 0.251, 0.213, 0.213, 0.214, 0.216, 0.216, 0.209,
      0.190, 0.181, 0.183, 0.187, 0.187, 0.183, 0.184, 0.193))
})

test_that("- pcvpcr() works well for Corn data with local scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # because of error in ordering of CV values (fixed now) we have to provide
   # manual vector in this test
   cv <- rep(seq_len(4), length.out = nrow(X))[order(Y)]

   Xpv <- pcvpcr(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "local")
   D <- attr(Xpv, "D")

   rcv <- pcrcvlocal(X, Y, ncomp = 20, cv = cv, center = TRUE, scale = FALSE)
   rpv <- pcrpv(X, Y, Xpv, ncomp = 20, center = TRUE, scale = FALSE)

   expect_equivalent(rcv$Yp, rpv$Yp, tolerance = 0.02)
   expect_equivalent(rcv$RMSE, rpv$RMSE, tolerance = 0.02)
})

test_that("- pcvpcr() compare local and global CV scope.", {
   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   cv <- list("ven", 10)
   Xpv1 <- pcvpcr(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "global")
   Xpv2 <- pcvpcr(X, Y, 20, center = TRUE, scale = FALSE, cv = cv, cv.scope = "local")

   rcv1 <- pcrcvglobal(X, Y, ncomp = 20, cv = cv, center = TRUE, scale = FALSE)
   rcv2 <- pcrcvlocal(X, Y, ncomp = 20, cv = cv, center = TRUE, scale = FALSE)
   rpv1 <- pcrpv(X, Y, Xpv1, ncomp = 20, center = TRUE, scale = FALSE)
   rpv2 <- pcrpv(X, Y, Xpv2, ncomp = 20, center = TRUE, scale = FALSE)

   expect_equivalent(rcv1$RMSE, rpv1$RMSE)
   expect_equivalent(rcv2$RMSE, rpv2$RMSE, tolerance = 0.01)

})
