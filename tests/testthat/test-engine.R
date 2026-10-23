small_problem <- function(n = 60, G = 8, seed = 1) {
  set.seed(seed)
  sim <- simulate_signal_data(n_cells = n, n_genes = c(null = G / 2, ctSVG = G / 2),
                              seed = seed)
  Y <- t(apply(sim$counts, 1, rank_int_impute))
  K_ct <- build_celltype_kernel(sim$celltype)
  K_sp <- build_spatial_kernels(sim$location)
  list(Y = Y, K_ct = K_ct, K_sp = K_sp, sim = sim)
}

test_that(".liu matches CompQuadForm::liu", {
  skip_if_not_installed("CompQuadForm")
  set.seed(2)
  for (i in 1:20) {
    lambda <- rnorm(30)
    q <- rnorm(1)
    expect_equal(SIGNAL:::.liu_eigen(q, lambda), CompQuadForm::liu(q, lambda))
  }
  lambda <- rexp(30)
  expect_equal(SIGNAL:::.liu_eigen(5, lambda), CompQuadForm::liu(5, lambda))
})

test_that("saddlepoint matches exact and Imhof tail probabilities", {
  sp <- SIGNAL:::.saddle_eigen
  for (q in c(5, 20, 80, 160)) {
    expect_equal(sp(q, rep(1, 10)), pchisq(q, 10, lower.tail = FALSE), tolerance = 0.02)
  }
  set.seed(1)
  lambda <- c(rexp(40), -rexp(20, 2))
  for (q in c(0, 30, 60, 100)) {
    expect_equal(sp(q, lambda), CompQuadForm::imhof(q, lambda, epsabs = 1e-14, epsrel = 1e-14,
                                                   limit = 1e4)$Qq, tolerance = 0.03)
  }
  expect_equal(sp(1, -rexp(5)), 0)
  expect_equal(sp(-1, rexp(5)), 1)
})

test_that("saddlepoint method agrees with davies on simulated genes", {
  pb <- small_problem()
  Klist <- SIGNAL:::.signal_kernels(pb$K_ct, pb$K_sp[[3]])
  p_dav <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "davies")
  p_sad <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "saddlepoint")
  ok <- p_dav > 1e-4
  expect_lt(max(abs(log10(p_sad[ok] / p_dav[ok]))), 0.1)
  expect_true(all(p_sad > 0))
})

test_that("kernel list matches the original construction", {
  pb <- small_problem()
  for (Ks in pb$K_sp[c(1, 6, 11)]) {
    expect_equal(SIGNAL:::.signal_kernels(pb$K_ct, Ks), ref_kernels(pb$K_ct, Ks))
  }
})

test_that("eigen method reproduces the original implementation exactly", {
  skip_if_not_installed("CompQuadForm")
  pb <- small_problem()
  for (Ks in pb$K_sp[c(1, 4, 11)]) {
    Klist <- SIGNAL:::.signal_kernels(pb$K_ct, Ks)
    ours <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "eigen")
    ref <- t(apply(pb$Y, 1, ref_test_gene, KList = Klist))
    expect_equal(ours, ref, ignore_attr = TRUE, tolerance = 1e-8)
  }
})

test_that("davies method reproduces the original tests with CompQuadForm::davies", {
  pb <- small_problem()
  davies_p <- function(q, lambda) min(max(CompQuadForm::davies(q, lambda)$Qq, 0), 1)
  for (Ks in pb$K_sp[c(1, 11)]) {
    Klist <- SIGNAL:::.signal_kernels(pb$K_ct, Ks)
    ours <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "davies")
    ref <- t(apply(pb$Y, 1, ref_test_gene, KList = Klist, pfun = davies_p))
    expect_equal(ours, ref, ignore_attr = TRUE, tolerance = 1e-8)
  }
})

test_that("davies_liu method reproduces the original tests (Davies overall, Liu components)", {
  pb <- small_problem()
  Yr <- t(apply(pb$sim$counts, 1, function(y) qnorm((rank(y) - 0.5) / length(y))))
  for (Ks in pb$K_sp[c(1, 6, 11)]) {
    Klist <- SIGNAL:::.signal_kernels(pb$K_ct, Ks)
    ours <- SIGNAL:::.signal_engine(Yr, Klist, 1:3, method = "davies_liu")
    ref <- t(apply(Yr, 1, ref_test_gene, KList = Klist, pfun = CompQuadForm::liu,
                   pfun_overall = function(q, lambda) CompQuadForm::davies(q, lambda)$Qq))
    expect_equal(ours, ref, ignore_attr = TRUE, tolerance = 1e-6)
  }
})

test_that("SIGNAL() with davies_liu reproduces the original real-data pipeline end to end", {
  sim <- simulate_signal_data(n_cells = 70, n_genes = c(null = 3, spatial = 2, ctSVG = 3), seed = 11)
  K <- build_spatial_kernels(sim$location)[c(1, 4, 11)]
  K_ct <- build_celltype_kernel(sim$celltype)
  # Original: plain ranks, Testing_Each_Gene per kernel, CombinePValues over kernels
  Yr <- t(apply(sim$counts, 1, function(y) qnorm((rank(y) - 0.5) / length(y))))
  per_kernel <- lapply(K, function(Ks) t(apply(Yr, 1, ref_test_gene,
    KList = ref_kernels(K_ct, Ks), pfun = CompQuadForm::liu,
    pfun_overall = function(q, lambda) CompQuadForm::davies(q, lambda)$Qq)))
  ref <- sapply(1:4, function(t) CombinePValues(sapply(per_kernel, function(p) p[, t])))
  obj <- SIGNAL(sim$counts, sim$location, sim$celltype, kernels = K, min_cells = 0, transform = "rank", method = "davies_liu",
                verbose = FALSE)
  ours <- as.matrix(obj@result[, c("p_overall", "p_celltype", "p_spatial", "p_ctSVG")])
  expect_equal(ours, ref, ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("liu method agrees with Liu on per-gene eigenvalues", {
  pb <- small_problem()
  for (Ks in pb$K_sp) {
    Klist <- SIGNAL:::.signal_kernels(pb$K_ct, Ks)
    p_eig <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "eigen")
    p_mom <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "liu")
    expect_equal(p_mom, p_eig, tolerance = 1e-6)
  }
})

test_that("spatial kernels follow their closed-form definitions", {
  set.seed(3)
  loc <- cbind(runif(80, 0, 10), runif(80, 0, 5))
  K <- build_spatial_kernels(loc)
  expect_named(K, c(paste0(rep(c("gaussian_", "matern_"), 5), rep(1:5, each = 2)),
                    "spline"))

  z <- scale(loc)
  d <- unname(as.matrix(dist(z)))
  l <- quantile(abs(z), probs = seq(0.2, 1, by = 0.2))
  for (s in 1:5) {
    expect_equal(K[[paste0("gaussian_", s)]], exp(-d^2 / (2 * l[[s]]^2)))
    expect_equal(K[[paste0("matern_", s)]], (1 + l[[s]] * d) * exp(-l[[s]] * d))
  }

  S <- K$spline
  expect_equal(S, t(S))
  expect_equal(S %*% S, S)
  expect_equal(sum(diag(S)), 4)
})

test_that("cell-type kernel uses Hamming for single-cell labels", {
  ct <- c("a", "b", "a", "c")
  K <- build_celltype_kernel(ct)
  expect_equal(diag(K), rep(1, 4))
  expect_equal(K[1, 3], 1)
  expect_equal(K[1, 2], exp(-2))
  # Multi-label 0/1 rows: Hamming counts differing columns (original hamming_kernel)
  X <- rbind(c(1, 0, 0), c(1, 1, 0), c(0, 1, 1))
  ref <- exp(-unname(as.matrix(stats::dist(X, method = "manhattan"))))
  expect_equal(build_celltype_kernel(X), ref)
  expect_equal(build_celltype_kernel(X, distance = "hamming"), ref)
})

test_that("cell-type kernel on proportions is the original Gaussian_kernel", {
  set.seed(8)
  P <- matrix(rexp(50 * 4), 50)
  P <- P / rowSums(P)
  # Gaussian_kernel() in funcs/ctfuncs.R of the original analysis code
  ref <- exp(-1 * as.matrix(stats::dist(P)^2))
  expect_equal(build_celltype_kernel(P), unname(ref))
  expect_equal(build_celltype_kernel(P, distance = "euclidean", gamma = 2), unname(ref^2))
  sim <- simulate_signal_data(n_cells = 50, n_genes = c(null = 4, spatial = 4), seed = 2)
  obj <- createSIGNALobject(sim$counts, sim$location, P)
  obj <- SIGNAL_kernel(SIGNAL_preprocess(obj, min_cells = 0, seed = 1, verbose = FALSE),
                       verbose = FALSE)
  expect_equal(obj@params$kernel$celltype_distance, "euclidean")
  expect_equal(obj@celltype_kernel, unname(ref))
})

test_that("combine_pvalues is the original CombinePValues", {
  set.seed(5)
  P <- matrix(runif(11 * 200), ncol = 11)
  P[1, 3] <- 0
  P[2, 5] <- 0.9999
  P[3, ] <- 1e-20
  expect_identical(combine_pvalues(P), CombinePValues(P))
  W <- matrix(c(rep(1 / 15, 10), 1 / 3), nrow = nrow(P), ncol = 11, byrow = TRUE)
  expect_identical(combine_pvalues(P, W), CombinePValues(P, W))
  P[4, 2] <- NA
  expect_identical(combine_pvalues(P), CombinePValues(P))
  expect_true(is.na(combine_pvalues(P)[4]))
  expect_equal(combine_pvalues(P, na.rm = TRUE)[4], combine_pvalues(P[4, -2, drop = FALSE]))
})

test_that("single Gaussian kernel follows CELINA's approximation kernel", {
  sim <- simulate_signal_data(n_cells = 80, n_genes = c(null = 5, spatial = 5), seed = 4)
  lib <- colSums(sim$counts)
  norm <- log2(t(t(sim$counts) / (lib / mean(lib))) + 1)
  bw <- median(apply(norm[1:6, ], 1, stats::bw.nrd0))
  expect_equal(select_bandwidth(sim$counts, genes = 1:6), bw)
  K <- build_gaussian_kernel(sim$location, counts = sim$counts, genes = 1:6, tol = 1e-4)
  ref <- exp(-unname(as.matrix(dist(scale(sim$location))))^2 / bw)
  ref[ref < 1e-4] <- 0
  expect_equal(c(K), c(ref))
  expect_equal(attr(K, "bandwidth"), bw)
})

test_that("SIGNAL_kernel switches to one kernel above max_cells", {
  sim <- simulate_signal_data(n_cells = 60, n_genes = c(null = 4, spatial = 4), seed = 6)
  obj <- SIGNAL_preprocess(createSIGNALobject(sim$counts, sim$location, sim$celltype),
                           min_cells = 0, seed = 1, verbose = FALSE)
  expect_length(SIGNAL_kernel(obj, verbose = FALSE)@spatial_kernels, 11)
  one <- SIGNAL_kernel(obj, max_cells = 50, verbose = FALSE)
  expect_named(one@spatial_kernels, "gaussian")
  expect_equal(one@params$kernel$bandwidth,
               select_bandwidth(sim$counts, rownames(obj@gene_expression)))
  expect_length(SIGNAL_kernel(obj, approximation = TRUE, verbose = FALSE)@spatial_kernels, 1)
  expect_length(SIGNAL_kernel(obj, approximation = FALSE, max_cells = 50,
                              verbose = FALSE)@spatial_kernels, 11)
  res <- SIGNAL_test(one, verbose = FALSE)@result
  expect_true(all(res$p_overall >= 0 & res$p_overall <= 1))
})
