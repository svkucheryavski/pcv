####################################
# Tests for PCV for three-way data #
####################################

# simulated three-way data: three trilinear components plus noise
simdata3w <- function(I = 24, J = 10, K = 8, seed = 1) {
   set.seed(seed)
   A <- matrix(rnorm(I * 3), I)
   B <- matrix(runif(J * 3), J)
   C <- matrix(runif(K * 3), K)
   X <- array(0, c(I, J, K))
   for (f in 1:3) X <- X + outer(outer(A[, f], B[, f]), C[, f])
   X + array(rnorm(I * J * K, sd = 0.05), dim(X)) + 5
}

# reference implementation, written independently of the package code
ref.prep <- function(X, center, scale) {
   M <- if (center) apply(X, c(2, 3), mean) else matrix(0, dim(X)[2], dim(X)[3])
   Xp <- sweep(X, c(2, 3), M, "-")
   S2 <- if (scale %in% c("mode2", "both")) apply(Xp, 2, sd) else rep(1, dim(X)[2])
   Xp <- sweep(Xp, 2, S2, "/")
   S3 <- if (scale %in% c("mode3", "both")) apply(Xp, 3, sd) else rep(1, dim(X)[3])
   Xp <- sweep(Xp, 3, S3, "/")
   list(Xp = Xp, apply = function(Y) sweep(sweep(sweep(Y, c(2, 3), M, "-"), 2, S2, "/"), 3, S3, "/"))
}

ref.factors <- function(X, mode, R) {
   d <- dim(X)
   M <- if (mode == 2) matrix(aperm(X, c(2, 3, 1)), d[2]) else matrix(aperm(X, c(3, 2, 1)), d[3])
   svd(M, nu = R, nv = 0)$u
}

ref.slab <- function(X, i) matrix(X[i, , ], dim(X)[2], dim(X)[3])

# check that projections of HSV PV-set on the global model equal cross-validation results
expect_hsv_cv <- function(X, P, center, scale, tol = 1e-8) {
   R <- attr(P, "ncomp")
   cvind <- attr(P, "cvind")
   al <- attr(P, "alignment")
   vs <- attr(P, "valid.subranks")

   p <- ref.prep(X, center, scale)
   Pp <- p$apply(P)
   U2 <- ref.factors(p$Xp, 2, R[1])
   U3 <- ref.factors(p$Xp, 3, R[2])

   for (k in unique(cvind)) {
      L2 <- ref.factors(p$Xp[cvind != k, , , drop = FALSE], 2, R[1])
      L3 <- ref.factors(p$Xp[cvind != k, , , drop = FALSE], 3, R[2])

      # local factors in the order reported by the method, signs matched to the global ones
      U2k <- L2[, al$local.col[al$segment == k & al$mode == 2], drop = FALSE]
      U3k <- L3[, al$local.col[al$segment == k & al$mode == 3], drop = FALSE]
      U2k <- U2k %*% diag(sign(colSums(U2 * U2k)), R[1], R[1])
      U3k <- U3k %*% diag(sign(colSums(U3 * U3k)), R[2], R[2])

      for (i in which(cvind == k)) {
         Xi <- ref.slab(p$Xp, i)
         Pi <- ref.slab(Pp, i)

         # all components: scores and orthogonal distance
         Tcv <- crossprod(U2k, Xi) %*% U3k
         Tpv <- crossprod(U2, Pi) %*% U3
         expect_equal(Tpv, Tcv, tolerance = tol)
         expect_equal(sum(Pi^2) - sum(Tpv^2), sum((Xi - U2k %*% Tcv %*% t(U3k))^2), tolerance = tol)

         # valid smaller numbers of components: compare with local models fitted directly
         for (r2 in vs$r2) for (r3 in vs$r3) {
            l2 <- L2[, seq_len(r2), drop = FALSE]
            l3 <- L3[, seq_len(r3), drop = FALSE]
            qcv <- sum((Xi - l2 %*% crossprod(l2, Xi) %*% l3 %*% t(l3))^2)
            qpv <- sum(Pi^2) - sum(Tpv[seq_len(r2), seq_len(r3)]^2)
            expect_equal(qpv, qcv, tolerance = tol)
         }
      }
   }
}


context("pcv3w: HSV method")

test_that("HSV PV-set reproduces cross-validation for all preprocessing settings", {
   X <- simdata3w()
   for (center in c(TRUE, FALSE)) {
      for (scale in c("none", "mode2", "mode3", "both")) {
         P <- suppressWarnings(pcv3whsv(X, c(3, 3), cv = list("ven", 4), center = center,
            scale = scale, seed = 1))
         expect_equal(dim(P), dim(X))
         expect_equal(attr(P, "ncomp"), c(3, 3))
         expect_hsv_cv(X, P, center, scale)
      }
   }
})

test_that("HSV PV-set reproduces cross-validation for different splits and settings", {
   X <- simdata3w(seed = 2)
   for (cv in list(list("ven", 3), list("rand", 5), list("loo"), rep(1:4, length.out = 24))) {
      for (orth in c("random", "canonical")) {
         for (align in c("optimal", "greedy")) {
            P <- suppressWarnings(pcv3whsv(X, c(4, 3), cv = cv, orth = orth, align = align, seed = 3))
            expect_hsv_cv(X, P, TRUE, "none")
         }
      }
   }
})

test_that("HSV PV-set is valid for smaller numbers of components", {
   X <- simdata3w(seed = 3)
   # catch the warning about invalid smaller numbers of components
   warned <- FALSE
   P <- withCallingHandlers(pcv3whsv(X, c(5, 5), cv = list("ven", 4), seed = 1),
      warning = function(w) {
         warned <<- grepl("smaller numbers of components", conditionMessage(w))
         invokeRestart("muffleWarning")
      })
   vs <- attr(P, "valid.subranks")

   # generation ranks are always valid, the warning comes only when some smaller are not
   expect_true(5 %in% vs$r2 && 5 %in% vs$r3)
   expect_equal(warned, length(vs$r2) < 5 || length(vs$r3) < 5)
   expect_hsv_cv(X, P, TRUE, "none")
})

test_that("HSV attributes are correct", {
   X <- simdata3w()
   P <- pcv3whsv(X, c(3, 2), cv = list("ven", 4), scale = TRUE, orth = "canonical", seed = 1)
   expect_equal(attr(P, "cvind"), rep(1:4, length.out = 24))
   expect_equal(attr(P, "orth"), "canonical")
   expect_equal(attr(P, "align"), "optimal")
   expect_equal(attr(P, "scale"), "both")
   expect_equal(dim(attr(P, "prep")$M), c(10, 8))
   expect_equal(nrow(attr(P, "alignment")), 4 * (3 + 2))
   expect_true(all(attr(P, "alignment")$cosine <= 1 + 1e-12))
})

test_that("HSV with canonical orthogonal part does not depend on random numbers", {
   X <- simdata3w()
   P1 <- pcv3whsv(X, c(3, 3), orth = "canonical", seed = 1)
   P2 <- pcv3whsv(X, c(3, 3), orth = "canonical", seed = 2)
   expect_equal(P1, P2)
})

test_that("HSV reduces too large number of components with a warning", {
   X <- simdata3w(I = 12, J = 4, K = 5)
   expect_warning(P <- pcv3whsv(X, c(6, 6), cv = list("ven", 3), seed = 1), "reduced")
   expect_equal(attr(P, "ncomp"), c(4, 5))
   expect_hsv_cv(X, P, TRUE, "none")
})


context("pcv3w: FU and BU methods")

test_that("FU PV-set is the same as pcvpca() on unfolded data", {
   X <- simdata3w()
   d <- dim(X)
   for (scale in c(FALSE, TRUE)) {
      P <- pcv3wfu(X, ncomp = 6, cv = list("ven", 4), scale = scale, seed = 5)
      set.seed(5)
      Xpv <- pcvpca(matrix(X, d[1]), ncomp = 6, cv = list("ven", 4), scale = scale)
      expect_equal(P, array(Xpv, d))
   }
})

test_that("BU PV-set is the same as pcvpca() on each block with one split", {
   X <- simdata3w()
   d <- dim(X)
   for (mode in 2:3) {
      P <- pcv3wbu(X, mode = mode, ncomp = 4, cv = list("rand", 4), seed = 5)
      set.seed(5)
      cvind <- pcvcrossval(list("rand", 4), d[1])[, 1]
      expect_equal(attr(P, "cvind"), cvind)
      for (b in seq_len(d[mode])) {
         Xb <- if (mode == 2) X[, b, ] else X[, , b]
         Pb <- if (mode == 2) P[, b, ] else P[, , b]
         expect_equal(Pb, pcvpca(Xb, ncomp = 4, cv = cvind))
      }
   }
})


context("pcv3w: seed and input checks")

test_that("seed makes results reproducible and keeps caller's random state", {
   X <- simdata3w()
   for (f in list(
      function(s) pcv3whsv(X, c(3, 3), seed = s),
      function(s) pcv3wfu(X, 5, seed = s),
      function(s) pcv3wbu(X, 2, 4, seed = s))) {
      set.seed(99)
      r1 <- runif(1)
      set.seed(99)
      P1 <- f(10)
      r2 <- runif(1)
      expect_equal(r1, r2)
      expect_equal(P1, f(10))
      expect_false(isTRUE(all.equal(P1, f(11))))
   }
})

test_that("wrong input gives an error", {
   X <- simdata3w()
   expect_error(pcv3whsv(matrix(1:20, 4), c(2, 2)), "three-way")
   expect_error(pcv3whsv(X[, , 1, drop = FALSE], c(2, 2)), "at least 2")
   expect_error(pcv3whsv(X, 2), "two positive integers")
   expect_error(pcv3whsv(X, c(0, 2)), "two positive integers")
   expect_error(pcv3whsv(X, c(2.5, 2)), "two positive integers")

   Xna <- X
   Xna[1, 1, 1] <- NA
   expect_error(pcv3whsv(Xna, c(2, 2)), "NA")

   Xc <- X
   Xc[, 2, ] <- 1
   expect_error(pcv3whsv(Xc, c(2, 2), scale = "mode2"), "zero variance")
   expect_error(pcv3whsv(X, c(2, 2), scale = "wrong"))

   expect_error(pcv3whsv(X, c(2, 2), cv.scope = "local"))
   expect_error(pcv3whsv(X, c(2, 2), press = TRUE))
   expect_error(pcv3wbu(X, mode = 1), "mode")
   expect_error(pcv3wfu(matrix(1:20, 4)), "three-way")
})


context("pcv3w: private test on Juices data")

# Juices data are not public: set PCV_JUICES to the path of juices.csv to run this test
test_that("HSV, FU and BU work on Juices data", {
   path <- Sys.getenv("PCV_JUICES")
   skip_if(path == "" || !file.exists(path), "Juices data not available (PCV_JUICES not set)")

   d <- read.csv(path, row.names = 1)
   X <- as.matrix(d[, -1])
   X3w <- array(0, c(nrow(X), 61, 81))
   for (i in seq_len(nrow(X))) X3w[i, , ] <- matrix(X[i, ], 61, 81)

   P <- pcv3whsv(X3w, c(10, 10), cv = list("ven", 3), seed = 1)
   expect_equal(attr(P, "ncomp"), c(10, 10))
   expect_hsv_cv(X3w, P, TRUE, "none", tol = 1e-6)

   Pfu <- pcv3wfu(X3w, ncomp = 20, cv = list("ven", 3), seed = 1)
   expect_equal(dim(Pfu), dim(X3w))
   expect_true(all(is.finite(Pfu)))

   for (mode in 2:3) {
      Pbu <- pcv3wbu(X3w, mode = mode, ncomp = 20, cv = list("ven", 3), seed = 1)
      expect_equal(dim(Pbu), dim(X3w))
      expect_true(all(is.finite(Pbu)))
   }
})
