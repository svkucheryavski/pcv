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

# all permutations of 1..n
ref.perms <- function(n) if (n == 1) list(1L) else
   unlist(lapply(ref.perms(n - 1), function(p) lapply(0:(n - 1), function(k) append(p, n, k))), FALSE)

# TRUE if assignment p gives the maximal sum of M[i, p[i]] (TRUE without checking above 6 x 6)
ref.isoptimal <- function(M, p) {
   n <- nrow(M)
   if (n > 6) return(TRUE)
   best <- max(sapply(ref.perms(n), function(pp) sum(M[cbind(seq_len(n), pp)])))
   abs(sum(M[cbind(seq_len(n), p)]) - best) < 1e-12
}

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

      # reported matching must reach the maximal sum of absolute cosines (checked up to 6 x 6)
      p2 <- al$local.col[al$segment == k & al$mode == 2]
      p3 <- al$local.col[al$segment == k & al$mode == 3]
      if (attr(P, "align") == "optimal") {
         expect_true(ref.isoptimal(abs(crossprod(U2, L2)), p2))
         expect_true(ref.isoptimal(abs(crossprod(U3, L3)), p3))
      }

      # local factors in the order reported by the method, signs matched to the global ones
      U2k <- L2[, p2, drop = FALSE]
      U3k <- L3[, p3, drop = FALSE]
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
   expect_true(length(vs$r2) < 5)
   expect_true(warned)
   expect_hsv_cv(X, P, TRUE, "none")

   # omitted numbers of components really do not reproduce cross-validation
   cvind <- attr(P, "cvind")
   Xc <- sweep(X, c(2, 3), apply(X, c(2, 3), mean), "-")
   Pc <- sweep(P, c(2, 3), apply(X, c(2, 3), mean), "-")
   U2 <- ref.factors(Xc, 2, 5)
   U3 <- ref.factors(Xc, 3, 5)
   for (mode in 2:3) {
      for (r in setdiff(1:5, vs[[mode - 1]])) {
         r2 <- if (mode == 2) r else 5
         r3 <- if (mode == 3) r else 5
         maxdiff <- 0
         for (k in unique(cvind)) {
            l2 <- ref.factors(Xc[cvind != k, , , drop = FALSE], 2, r2)
            l3 <- ref.factors(Xc[cvind != k, , , drop = FALSE], 3, r3)
            for (i in which(cvind == k)) {
               Xi <- ref.slab(Xc, i)
               Pi <- ref.slab(Pc, i)
               qcv <- sum((Xi - l2 %*% crossprod(l2, Xi) %*% l3 %*% t(l3))^2)
               qpv <- sum(Pi^2) - sum((crossprod(U2[, seq_len(r2), drop = FALSE], Pi) %*%
                  U3[, seq_len(r3), drop = FALSE])^2)
               maxdiff <- max(maxdiff, abs(qpv - qcv) / qcv)
            }
         }
         expect_gt(maxdiff, 1e-6)
      }
   }
})

test_that("optimal and greedy matching give expected assignments", {
   # absolute cosines (global components in rows, local in columns) where the methods differ
   M <- matrix(c(
      0.90, 0.80, 0.10,
      0.85, 0.10, 0.30,
      0.20, 0.75, 0.70), 3, 3, byrow = TRUE)
   expect_equal(lsap3w(M), c(2L, 1L, 3L))

   Ug <- diag(3)
   Ul <- M * matrix(c(1, -1, 1), 3, 3, byrow = TRUE)
   expect_equal(attr(align3w(Ug, Ul, "optimal"), "perm"), c(2L, 1L, 3L))
   expect_equal(attr(align3w(Ug, Ul, "greedy"), "perm"), c(1L, 3L, 2L))

   # signs are flipped so matched pairs have positive inner products
   Ua <- align3w(Ug, Ul, "optimal")
   expect_true(all(diag(crossprod(Ug, Ua)) > 0))
   expect_equal(attr(Ua, "cos"), c(0.80, 0.85, 0.70))

   # optimal matching finds the maximal sum for random matrices up to 6 x 6
   perms <- function(n) if (n == 1) list(1L) else
      unlist(lapply(perms(n - 1), function(p) lapply(0:(n - 1), function(k) append(p, n, k))), FALSE)
   set.seed(7)
   for (it in 1:100) {
      n <- sample(2:6, 1)
      M <- matrix(runif(n * n), n)
      best <- max(sapply(perms(n), function(p) sum(M[cbind(seq_len(n), p)])))
      expect_equal(sum(M[cbind(seq_len(n), lsap3w(M))]), best)
   }
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

      # column (j, k) of the unfolded data is column j + J * (k - 1)
      for (j in c(1, d[2])) for (k in c(1, 2, d[3])) {
         expect_equal(P[, j, k], Xpv[, j + d[2] * (k - 1)])
      }
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


test_that("wrong segment numbers in 'cv' give an error", {
   X <- simdata3w()
   for (cv in list(rep(0:2, 8), rep(c(1, 3), 12), rep(1, 24), rep(c(1.5, 2), 12))) {
      expect_error(pcv3whsv(X, c(2, 2), cv = cv), "segment numbers")
      expect_error(pcv3wbu(X, 2, 3, cv = cv), "segment numbers")
      expect_error(pcv3wfu(X, 3, cv = cv), "segment numbers")
   }
   expect_warning(pcv3wfu(X, 3, cv = list("rand", 4, 3), seed = 1), "first one")
})

test_that("orthogonal part uses the fixed direction when needed", {
   set.seed(5)
   U2 <- qr.Q(qr(matrix(rnorm(6 * 2), 6)))
   U3 <- qr.Q(qr(matrix(rnorm(5 * 2), 5)))

   # direction does not depend on signs of the factors
   cdir <- complementdir3w(U2, U3)
   expect_equal(complementdir3w(-U2, -U3), cdir)
   expect_equal(sum(cdir^2), 1)
   expect_equal(sum(proj3w(cdir, U2, U3)^2), 0)

   # residual lying inside the model space: the fixed direction is used with the right norm
   E <- U2 %*% matrix(c(1, 2, 3, 4), 2) %*% t(U3)
   for (orth in c("canonical", "random")) {
      Xk <- if (orth == "canonical") list(E) else list(E, 2 * E)
      Eres <- if (orth == "canonical") list(E) else list(E, E)
      Eo <- getxpvorth3w(Xk, Eres, U2, U3, cdir, orth)
      expect_equal(Eo[[1]], cdir * sqrt(sum(E^2)))
   }

   # full model space in both modes: no direction
   expect_null(complementdir3w(diag(3), diag(2)))
})

test_that("names of dimensions are kept", {
   X <- simdata3w()
   dimnames(X) <- list(paste0("s", 1:24), paste0("ex", 1:10), paste0("em", 1:8))
   expect_equal(dimnames(pcv3whsv(X, c(2, 2), seed = 1)), dimnames(X))
   expect_equal(dimnames(pcv3wfu(X, 3, seed = 1)), dimnames(X))
   expect_equal(dimnames(pcv3wbu(X, 2, 3, seed = 1)), dimnames(X))
})


context("pcv3w: reference values for Python and MATLAB")

## directories with dataset and reference values shared with Python and MATLAB tests
dataFile3w <- "../../../.tests/data/sim3w-20x8x6.csv"
caseDir3w <- "../../../.tests/pcv3whsv/"

# simulated data for references: Tucker structure with random 4 x 3 core per sample plus noise
refdata3w <- function() {
   set.seed(42)
   I <- 20; J <- 8; K <- 6
   B <- qr.Q(qr(matrix(rnorm(J * 4), J)))
   C <- qr.Q(qr(matrix(rnorm(K * 3), K)))
   X <- array(0, c(I, J, K))
   for (i in seq_len(I)) X[i, , ] <- B %*% matrix(rnorm(12, sd = c(4, 3, 2, 1)), 4, 3) %*% t(C)
   round(X + array(rnorm(I * J * K, sd = 0.1), dim(X)) + 3, 6)
}

# q and h of a PV-set on the global model (h uses the pseudo-inverse of the score scatter)
ref.qh <- function(X, P, R, center, scale) {
   p <- ref.prep(X, center, scale)
   U2 <- ref.factors(p$Xp, 2, R[1])
   U3 <- ref.factors(p$Xp, 3, R[2])
   scores <- function(A) t(sapply(seq_len(dim(A)[1]), function(i) as.vector(crossprod(U2, ref.slab(A, i)) %*% U3)))
   s <- svd(scores(p$Xp), nu = 0)
   use <- s$d > sqrt(.Machine$double.eps) * s$d[1]
   Pp <- p$apply(P)
   T <- scores(Pp)
   q <- sapply(seq_len(dim(P)[1]), function(i) sum(ref.slab(Pp, i)^2)) - rowSums(T^2)
   h <- rowSums((T %*% s$v[, use, drop = FALSE] %*% diag(1 / s$d[use], sum(use)))^2)
   list(q = q, h = h)
}

# test cases: R2, R3, center, scale, align, cv
refcases3w <- function() {
   cases <- list()
   for (center in c(TRUE, FALSE)) for (scale in c("none", "mode2", "mode3", "both")) {
      cases[[length(cases) + 1]] <- list(ncomp = c(3, 3), center = center, scale = scale,
         align = "optimal", cv = list("ven", 4))
   }
   for (align in c("optimal", "greedy")) for (cv in list(list("ven", 4), list("loo"))) {
      cases[[length(cases) + 1]] <- list(ncomp = c(4, 2), center = TRUE, scale = "none",
         align = align, cv = cv)
   }
   cases
}

refsuffix3w <- function(cs) {
   cvtext <- if (length(cs$cv) == 2) paste0(cs$cv[[1]], cs$cv[[2]]) else cs$cv[[1]]
   sprintf("-%dx%d-%s-%s-%s-%s.csv", cs$ncomp[1], cs$ncomp[2], toupper(cs$center), cs$scale,
      cs$align, cvtext)
}

# set PCV_SAVE_REFS=TRUE to (re)generate the dataset and the reference values
test_that("HSV results match the reference values", {
   save <- Sys.getenv("PCV_SAVE_REFS") == "TRUE"
   wcsv <- function(x, f) write.table(x, f, col.names = FALSE, row.names = FALSE, sep = ",", dec = ".")
   rcsv <- function(f) unname(as.matrix(read.csv(f, header = FALSE)))

   if (save) {
      dir.create(caseDir3w, showWarnings = FALSE)
      X <- refdata3w()
      wcsv(matrix(X, dim(X)[1]), dataFile3w)
   }
   skip_if_not(file.exists(dataFile3w), "reference dataset not found")

   # dataset is stored unfolded, I x (J*K) with mode 2 running fastest
   X <- array(rcsv(dataFile3w), c(20, 8, 6))

   for (cs in refcases3w()) {
      sfx <- refsuffix3w(cs)

      # canonical orthogonal part: PV-set itself is fully determined
      P <- suppressWarnings(pcv3whsv(X, cs$ncomp, cv = cs$cv, center = cs$center, scale = cs$scale,
         orth = "canonical", align = cs$align))
      al <- unname(as.matrix(attr(P, "alignment")))

      # random orthogonal part: q and h are still determined (equal to cross-validation)
      Pr <- suppressWarnings(pcv3whsv(X, cs$ncomp, cv = cs$cv, center = cs$center, scale = cs$scale,
         orth = "random", align = cs$align, seed = 1))
      qh <- ref.qh(X, Pr, attr(Pr, "ncomp"), cs$center, cs$scale)
      qhc <- ref.qh(X, P, attr(P, "ncomp"), cs$center, cs$scale)
      expect_equal(qh$q, qhc$q, tolerance = 1e-8)
      expect_equal(qh$h, qhc$h, tolerance = 1e-8)

      if (save) {
         wcsv(matrix(P, dim(P)[1]), paste0(caseDir3w, "Xpv", sfx))
         wcsv(al, paste0(caseDir3w, "align", sfx))
         wcsv(qh$q, paste0(caseDir3w, "Q", sfx))
         wcsv(qh$h, paste0(caseDir3w, "H", sfx))
      }

      expect_equal(matrix(P, dim(P)[1]), rcsv(paste0(caseDir3w, "Xpv", sfx)), tolerance = 1e-8)
      expect_equal(al, rcsv(paste0(caseDir3w, "align", sfx)), tolerance = 1e-8)
      expect_equal(qh$q, as.vector(rcsv(paste0(caseDir3w, "Q", sfx))), tolerance = 1e-8)
      expect_equal(qh$h, as.vector(rcsv(paste0(caseDir3w, "H", sfx))), tolerance = 1e-8)
   }
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
