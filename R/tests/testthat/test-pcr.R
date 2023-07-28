####################################
# Tests for PCV for PCR            #
####################################

## directory to keep CSV files with reference values
caseDir <- "../../../.tests/pcvpcr/"

setup({
   pdf(file = tempfile("pcv-test-pcr-", fileext = ".pdf"))
})

teardown({
   dev.off()
})

#' Show distance plot for the results
plotpredictions <- function(Y, rcl, rcv, rpv, a, main = "Predictions") {

   y <- Y[, 1]
   yp <- rcl$Yp[, a]
   ypv <- rpv$Yp[, a]
   ycv <- rcv$Yp[, a]

   plot(
      c(y, y, y), c(yp, ycv, ypv),
      col = rep(c("blue", "red", "black"), each = length(y)),
      pch = rep(c(5, 1, 4), each = length(y)),
      xlab = "Reference values, y",
      ylab = "Predicted values, yp",
      main = main
   )

}

#' Create global PCR model and apply it to Xpv set
pcrpv <- function(X, Y, Xpv, ncomp, center = TRUE, scale = FALSE) {

   if (is.null(dim(Y))) {
      dim(Y) <- c(nrow(X), 1)
   }

   if (center) {
      mX <- apply(X, 2, mean)
      mY <- apply(Y, 2, mean)
   } else {
      mX <- rep(0, ncol(X))
      mY <- rep(0, ncol(Y))
   }

   if (scale)  {
      sX <- apply(X, 2, sd)
      sY <- apply(Y, 2, sd)
   } else {
      sX <- rep(1, ncol(X))
      sY <- rep(1, ncol(Y))
   }

   Xs <- scale(X, center = mX, scale = sX)
   Ys <- scale(Y, center = mY, scale = sY)
   Xpvs <- scale(Xpv, center = mX, scale = sX)

   # global model
   m <- svd(Xs, nu = ncomp, nv = ncomp)
   P <- m$v[, seq_len(ncomp), drop = FALSE]
   T <- Xs %*% P;

   # predictions for PV-set
   Tpv <- Xpvs %*% P;
   Ypv <- matrix(0, nrow(Y), ncomp)

   for (a in seq_len(ncomp)) {
      aind <- seq_len(a)
      Ta <- T[, aind, drop = FALSE]
      Tpva <- Tpv[, aind, drop = FALSE]
      Ca <- solve(crossprod(Ta)) %*% crossprod(Ta, Ys)
      Ypv[, a] <- Tpva %*% Ca;
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

   if (center) {
      mX <- apply(X, 2, mean)
      mY <- apply(Y, 2, mean)
   } else {
      mX <- rep(0, ncol(X))
      mY <- 0
   }

   if (scale) {
      sX <- apply(X, 2, sd)
      sY <- apply(Y, 2, sd)
   } else {
      sX <- rep(1, ncol(X))
      sY <- 1
   }

   Xs <- scale(X, center = mX, scale = sX)
   Ys <- scale(Y, center = mY, scale = sY)

   nSeg <- max(cvind)
   Ycv <- matrix(0, nrow(Y), ncomp)

   for (i in seq_len(nSeg)) {

      indc <- cvind != i
      indk <- cvind == i

      Xc <- Xs[indc, , drop = FALSE]
      Yc <- Ys[indc, , drop = FALSE]
      Xk <- Xs[indk, , drop = FALSE]

      mk <- svd(Xc, nu = ncomp, nv = ncomp)
      Pk <- mk$v[, seq_len(ncomp), drop = FALSE]

      Tc <- Xc %*% Pk;
      Tk <- Xk %*% Pk;

      for (a in seq_len(ncomp)) {
         aind <- seq_len(a)
         Tca <- Tc[, aind, drop = FALSE]
         Tka <- Tk[, aind, drop = FALSE]
         Ca <- solve(crossprod(Tca)) %*% crossprod(Tca, Yc)
         Ycv[indk, a] <- Tka %*% Ca;
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
   for (i in seq_len(nSeg)) {

      indc <- cvind != i
      indk <- cvind == i

      Xc <- X[indc, , drop = FALSE]
      Yc <- Y[indc, , drop = FALSE]
      Xk <- X[indk, , drop = FALSE]

      mX <- rep(0, ncol(Xc))
      mY <- 0
      if (center) {
         mX <- apply(Xc, 2, mean)
         mY <- apply(Yc, 2, mean)
      }

      sX <- rep(1, ncol(Xc))
      sY <- 1
      if (scale) {
         sX <- apply(Xc, 2, sd)
         sY <- apply(Yc, 2, sd)
      }


      Xcs <- scale(Xc, center = mX, scale = sX)
      Xks <- scale(Xk, center = mX, scale = sX)
      Ycs <- scale(Yc, center = mY, scale = sY)

      mk <- svd(Xcs, nu = ncomp, nv = ncomp)
      Pk <- mk$v[, seq_len(ncomp), drop = FALSE]
      Tc <- Xcs %*% Pk;
      Tk <- Xks %*% Pk;

      for (a in seq_len(ncomp)) {
         aind <- seq_len(a)
         Tca <- Tc[, aind, drop = FALSE]
         Tka <- Tk[, aind, drop = FALSE]
         Ca <- solve(crossprod(Tca)) %*% crossprod(Tca, Ycs)
         Ycv[indk, a] <- Tka %*% Ca;
      }

      Ycv[indk, ] <- Ycv[indk, ] * sY + mY;
   }

   E <- apply(Ycv, 2, function(v) v - Y)

   return (list(Yp = Ycv, RMSE = sqrt(colSums(E^2) / nrow(X))))
}

#' Save results as reference values
savereferences <- function(rpv.g, rpv.l, Dg, Dl, scale, ncomp, cv) {
   cvText <- if (length(cv) == 2) paste0(cv[[1]], cv[[2]]) else cv[[1]]
   caseSuffix <- sprintf("-%d-%s-%s.csv", ncomp, scale, cvText)

   write.table(rpv.g$Yp, file = paste0(caseDir, "Ypvg", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(rpv.l$Yp, file = paste0(caseDir, "Ypvl", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(Dg, file = paste0(caseDir, "Dg", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   write.table(Dl, file = paste0(caseDir, "Dl", caseSuffix), col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
}

#' Run tests for different combinations of parameters
runtests <- function(X, Y, tolerance = 10^-6, save.res = FALSE) {

   cv_cases <- list(list("ven", 4), list("ven", 10), list("loo"), list("rand", 10))
   ncomp_cases <- c(1, 10, 20, 30)
   scale_cases <- c(TRUE, FALSE)
   center <- TRUE

   cases <- expand.grid(cv = cv_cases, ncomp = ncomp_cases, scale = scale_cases)
   for (i in seq_len(nrow(cases))) {

      cv <- cases[i, "cv"][[1]]
      ncomp <- cases[i, "ncomp"]
      scale <- cases[i, "scale"]

      # test global results
      set.seed(42)
      Xpv.g <- pcvpcr(X, Y, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "global")

      set.seed(42)
      rcv.g <- pcrcvglobal(X, Y, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv.g <- pcrpv(X, Y, Xpv.g, ncomp = ncomp, center = center, scale = scale)

      # here we expect that the results are almost identical
      expect_equivalent(rcv.g$Yp, rpv.g$Yp, tolerance = tolerance)
      expect_equivalent(rcv.g$RMSE, rpv.g$RMSE, tolerance = tolerance)

      # test local results
      set.seed(42)
      Xpv.l <- pcvpcr(X, Y, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")

      set.seed(42)
      rcv.l <- pcrcvlocal(X, Y, ncomp = ncomp, cv = cv, center = center, scale = scale)
      rpv.l <- pcrpv(X, Y, Xpv.l, ncomp = ncomp, center = center, scale = scale)

      # compute intercept, slope and correlation between PV and CV predicted y-values
      err <- sapply(seq_len(ncomp), function(a) {
         m <- lm(rcv.l$Yp[, a] ~ rpv.l$Yp[, a])
         c(coefficients(m), cor(rcv.l$Yp[, a], rpv.l$Yp[, a]))
      })

      # here we expect a good match — slope = 1±0.03 and correlation >= 0.97
      expect_equivalent(abs(err[1, ]), rep(0, ncomp), tolerance = 0.50)
      expect_equivalent(err[2, ], rep(1, ncomp), tolerance = 0.03)
      expect_equivalent(err[3, ], rep(1, ncomp), tolerance = 0.03)

      # test scalars
      Dg <- attr(Xpv.g, "D")
      expect_false(is.null(Dg))
      expect_equal(ncol(Dg), ncomp)

      Dl <- attr(Xpv.l, "D")
      expect_false(is.null(Dl))
      expect_equal(ncol(Dl), ncomp)

      # save outcomes as reference
      if (save.res && cv[[1]] != "rand") {
         savereferences(rpv.g, rpv.l, Dg, Dl, scale, ncomp, cv)
      }

   }
}


context("Tests for 'pcvpcr()':")

test_that("automatic tests for Simdata data are passed.", {

   load("Simdata.RData")
   X <- simdata$spectra.c
   Y <- simdata$conc.c[, 1]

   runtests(X, Y, tolerance = 10^-4)
})

test_that("automatic tests for Corn data are passed.", {

   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # automatic tests
   runtests(X, Y, save.res = TRUE)
})

test_that("manual tests for Corn data are passed.", {

   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   # manual tests
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

test_that("visual tests for Corn data looks good.", {

   data(corn)
   X <- corn$spectra
   Y <- corn$moisture

   ncomp <- 30
   center <- TRUE
   scale <- TRUE
   cv <- list("ven", 4)

   # global scope
   Xpv.g <- pcvpcr(X, Y, ncomp = ncomp, center = center, scale = scale, cv = cv)
   rcv.g <- pcrcvglobal(X, Y, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rpv.g <- pcrpv(X, Y, Xpv.g, ncomp = ncomp, center = center, scale = scale)
   rcl.g <- pcrpv(X, Y, X, ncomp = ncomp, center = center, scale = scale)

   # local scope
   Xpv.l <- pcvpcr(X, Y, ncomp = ncomp, center = center, scale = scale, cv = cv, cv.scope = "local")
   rcv.l <- pcrcvlocal(X, Y, ncomp = ncomp, cv = cv, center = center, scale = scale)
   rpv.l <- pcrpv(X, Y, Xpv.l, ncomp = ncomp, center = center, scale = scale)
   rcl.l <- pcrpv(X, Y, X, ncomp = ncomp, center = center, scale = scale)

   # predictions plots
   par(mfrow = c(2, 2))
   plotpredictions(Y, rcl.g, rcv.g, rpv.g,  4, main = "Global, a = 4")
   plotpredictions(Y, rcl.g, rcv.g, rpv.g, 20, main = "Global, a = 20")
   plotpredictions(Y, rcl.l, rcv.l, rpv.l,  4, main = "Local, a = 4")
   plotpredictions(Y, rcl.l, rcv.l, rpv.l, 20, main = "Local, a = 20")

})
