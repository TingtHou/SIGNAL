test_that("SIGNAL returns a genes x p-values table that detects the example genes", {
  data(SIGNAL_example, package = "SIGNAL", envir = environment())
  obj <- SIGNAL(SIGNAL_example$counts, SIGNAL_example$location,
                SIGNAL_example$celltype, seed = 1, verbose = FALSE)
  expect_s4_class(obj, "SIGNAL")
  r <- obj@result
  expect_equal(nrow(r), 70)
  expect_named(r, c("p_overall", "p_celltype", "p_spatial", "p_ctSVG"))
  expect_true(all(as.matrix(r) >= 0 & as.matrix(r) <= 1))
  expect_equal(dim(obj@pvalues_by_kernel$ctSVG), c(70, 11))

  type <- SIGNAL_example$gene_type[rownames(r)]
  expect_gt(mean(r$p_overall[type != "null"] < 0.05), 0.8)
  expect_lt(mean(r$p_overall[type == "null"] < 0.05), 0.2)
  expect_gt(mean(r$p_celltype[type == "celltype"] < 0.05), 0.8)
  expect_gt(mean(r$p_spatial[type == "spatial"] < 0.05), 0.8)
  expect_gt(mean(r$p_ctSVG[type == "ctSVG"] < 0.05), 0.6)

  expect_output(show(obj), "SIGNAL object")
})

test_that("step-by-step workflow equals the wrapper and is reproducible", {
  sim <- simulate_signal_data(n_cells = 80, n_genes = c(null = 4, ctSVG = 4), seed = 5)
  K <- build_spatial_kernels(sim$location)[1:2]
  obj <- createSIGNALobject(sim$counts, sim$location, sim$celltype)
  obj <- SIGNAL_preprocess(obj, seed = 9, verbose = FALSE)
  obj <- SIGNAL_kernel(obj, kernels = K, verbose = FALSE)
  obj <- SIGNAL_test(obj, verbose = FALSE)
  wrap <- SIGNAL(sim$counts, sim$location, sim$celltype, kernels = K, seed = 9,
                 verbose = FALSE)
  expect_identical(obj@result, wrap@result)
  expect_equal(ncol(obj@pvalues_by_kernel$overall), 2)
})

test_that("SIGNAL_celltype finds the cell type of each ctSVG", {
  data(SIGNAL_example, package = "SIGNAL", envir = environment())
  obj <- createSIGNALobject(SIGNAL_example$counts, SIGNAL_example$location,
                            SIGNAL_example$celltype)
  obj <- SIGNAL_preprocess(obj, verbose = FALSE)
  obj <- SIGNAL_kernel(obj, verbose = FALSE)
  obj@spatial_kernels <- obj@spatial_kernels[c(1, 3)]
  genes <- c("gene1", "gene61", "gene62")

  # Single-cell level: the spatial test within each cell type
  sc <- SIGNAL_celltype(obj, genes = genes, ncores = 1, verbose = FALSE)
  expect_equal(sc@params$celltype$level, "single_cell")
  expect_equal(dim(sc@celltype_result), c(3, 4))
  expect_true(all(sc@celltype_result[c("gene61", "gene62"), "CT1"] < 1e-3))
  cells <- obj@celltype[, "CT1"] == 1
  y <- qnorm((rank(SIGNAL_example$counts["gene61", cells]) - 0.5) / sum(cells))
  Ks <- build_spatial_kernels(obj@location[cells, ])   # rebuilt for the subset
  p <- vapply(Ks, function(K) SIGNAL:::.signal_engine(t(y), list(K, diag(sum(cells))),
                                                       integer(0))[, 1], numeric(1))
  expect_equal(sc@celltype_result["gene61", "CT1"], combine_pvalues(t(p)))

  # Spot level: one interaction kernel per cell type
  sp <- SIGNAL_celltype(obj, genes = genes, level = "spot", ncores = 1, verbose = FALSE)
  expect_equal(dim(sp@celltype_result), c(3, 4))
  expect_true(all(sp@celltype_result[c("gene61", "gene62"), "CT1"] < 1e-3))
})

test_that("sparse counts give the same result as dense counts", {
  skip_if_not_installed("Matrix")
  sim <- simulate_signal_data(n_cells = 60, n_genes = c(null = 3), seed = 8)
  dense <- SIGNAL_preprocess(createSIGNALobject(sim$counts, sim$location, sim$celltype),
                             seed = 1, verbose = FALSE)
  sparse <- SIGNAL_preprocess(createSIGNALobject(Matrix::Matrix(sim$counts, sparse = TRUE),
                                                 sim$location, sim$celltype),
                              seed = 1, verbose = FALSE)
  expect_equal(sparse@gene_expression, dense@gene_expression)
  expect_equal(select_bandwidth(sparse@counts), select_bandwidth(dense@counts))
  expect_equal(select_bandwidth(sparse@counts, genes = 1:2), select_bandwidth(dense@counts, genes = 1:2))
})

test_that("filter = 'celina' matches CELINA's FilteringGenesCells", {
  FilteringGenesCells <- function(counts, percentage_cell = 0.1, min_total_counts = 10) {
    idx <- which(apply(counts, 1, function(x) sum(x != 0)) >= floor(percentage_cell * ncol(counts)))
    counts <- counts[idx, ]
    counts[, which(apply(counts, 2, sum) > min_total_counts)]
  }
  sim <- simulate_signal_data(n_cells = 120, n_genes = c(null = 30), seed = 3)
  counts <- sim$counts
  counts[1:5, ] <- 0
  counts[1:5, 1:5] <- 1                     # 5 genes expressed in only 5 cells
  counts[, 1:4] <- 0
  counts[6, 1:4] <- 1                       # 4 cells with total count 1
  ref <- FilteringGenesCells(counts, percentage_cell = 0.1, min_total_counts = 10)
  obj <- createSIGNALobject(counts, sim$location, sim$celltype)
  obj <- SIGNAL_preprocess(obj, filter = "celina", seed = 1, verbose = FALSE)
  expect_identical(rownames(obj@gene_expression), rownames(ref))
  expect_identical(colnames(obj@gene_expression), colnames(ref))
  expect_equal(nrow(obj@location), ncol(ref))
  expect_equal(nrow(obj@celltype), ncol(ref))
  expect_equal(ncol(obj@counts), ncol(ref))
})

test_that("input checks catch mismatched dimensions and missing steps", {
  sim <- simulate_signal_data(n_cells = 50, n_genes = c(null = 2), seed = 7)
  expect_error(createSIGNALobject(sim$counts, sim$location[-1, ], sim$celltype), "location")
  expect_error(createSIGNALobject(sim$counts, sim$location, sim$celltype[-1]), "celltype")
  obj <- createSIGNALobject(sim$counts, sim$location, sim$celltype)
  expect_error(SIGNAL_test(obj), "SIGNAL_preprocess")
})

test_that("SIGNAL_plotGene returns a ggplot", {
  data(SIGNAL_example, package = "SIGNAL", envir = environment())
  obj <- createSIGNALobject(SIGNAL_example$counts, SIGNAL_example$location,
                            SIGNAL_example$celltype)
  expect_s3_class(SIGNAL_plotGene(obj, "gene61", by_celltype = TRUE), "ggplot")
})

test_that("overall and component steps give the same p-values as SIGNAL_test", {
  sim <- simulate_signal_data(n_cells = 80, n_genes = c(null = 3, spatial = 3, ctSVG = 3), seed = 7)
  obj <- createSIGNALobject(sim$counts, sim$location, sim$celltype)
  obj <- SIGNAL_preprocess(obj, min_cells = 0, verbose = FALSE)
  obj <- SIGNAL_kernel(obj, kernels = build_spatial_kernels(sim$location)[1:2], verbose = FALSE)
  full <- SIGNAL_test(obj, ncores = 1, verbose = FALSE)

  step1 <- SIGNAL_overall(obj, ncores = 1, verbose = FALSE)
  expect_equal(step1@result$p_overall, full@result$p_overall)
  expect_true(all(is.na(as.matrix(step1@result[, c("p_celltype", "p_spatial", "p_ctSVG")]))))

  genes <- rownames(obj@gene_expression)[c(2, 5, 8)]
  step2 <- SIGNAL_component(step1, genes = genes, ncores = 1, verbose = FALSE)
  expect_equal(step2@result$p_overall, full@result$p_overall)
  expect_equal(step2@result[genes, ], full@result[genes, ])
  other <- setdiff(rownames(obj@gene_expression), genes)
  expect_true(all(is.na(as.matrix(step2@result[other, c("p_celltype", "p_spatial", "p_ctSVG")]))))
  expect_equal(step2@pvalues_by_kernel$ctSVG[genes, ], full@pvalues_by_kernel$ctSVG[genes, ])

  both <- SIGNAL_component(SIGNAL_overall(obj, ncores = 1, verbose = FALSE), ncores = 1, verbose = FALSE)
  expect_equal(both@result, full@result)
})
