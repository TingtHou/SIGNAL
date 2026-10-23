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

test_that("SIGNAL_celltype returns a genes x cell types matrix", {
  sim <- simulate_signal_data(n_cells = 80, n_genes = c(null = 2, ctSVG = 2), seed = 6)
  obj <- createSIGNALobject(sim$counts, sim$location, sim$celltype)
  obj <- SIGNAL_preprocess(obj, seed = 1, verbose = FALSE)
  obj <- SIGNAL_kernel(obj, kernels = build_spatial_kernels(sim$location)[1:2],
                       verbose = FALSE)
  obj <- SIGNAL_celltype(obj, genes = c("gene3", "gene4"), verbose = FALSE)
  expect_equal(dim(obj@celltype_result), c(2, 4))
  expect_equal(colnames(obj@celltype_result), levels(sim$celltype))
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
