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

test_that("SIGNAL() defaults reproduce the original pipeline end to end", {
  sim <- simulate_signal_data(n_cells = 70, n_genes = c(null = 3, spatial = 2, ctSVG = 3), seed = 11)
  K <- build_spatial_kernels(sim$location)[c(1, 4, 11)]
  K_ct <- build_celltype_kernel(sim$celltype)
  # Original: plain ranks, Testing_Each_Gene per kernel, CombinePValues over kernels
  Yr <- t(apply(sim$counts, 1, function(y) qnorm((rank(y) - 0.5) / length(y))))
  per_kernel <- lapply(K, function(Ks) t(apply(Yr, 1, ref_test_gene,
    KList = ref_kernels(K_ct, Ks), pfun = CompQuadForm::liu,
    pfun_overall = function(q, lambda) CompQuadForm::davies(q, lambda)$Qq)))
  ref <- sapply(1:4, function(t) CombinePValues(sapply(per_kernel, function(p) p[, t])))
  obj <- SIGNAL(sim$counts, sim$location, sim$celltype, kernels = K, min_cells = 0,
                verbose = FALSE)
  ours <- as.matrix(obj@result[, c("p_overall", "p_celltype", "p_spatial", "p_ctSVG")])
  expect_equal(ours, ref, ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("moment method agrees with the eigen method", {
  pb <- small_problem()
  for (Ks in pb$K_sp) {
    Klist <- SIGNAL:::.signal_kernels(pb$K_ct, Ks)
    p_eig <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "eigen")
    p_mom <- SIGNAL:::.signal_engine(pb$Y, Klist, 1:3, method = "moment")
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

test_that("cell-type kernel equals the Hamming kernel", {
  ct <- c("a", "b", "a", "c")
  K <- build_celltype_kernel(ct)
  expect_equal(diag(K), rep(1, 4))
  expect_equal(K[1, 3], 1)
  expect_equal(K[1, 2], exp(-2))
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
