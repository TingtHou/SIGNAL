#' Create a SIGNAL object
#'
#' @param counts A genes x cells matrix of raw counts (base matrix, data frame
#'   or sparse `dgCMatrix`). Row names are used as gene names and column names
#'   as cell names.
#' @param location A cells x 2 matrix or data frame of spatial coordinates.
#' @param celltype Cell-type labels (one per cell), or a cells x cell-types
#'   indicator matrix.
#' @param project Project name.
#' @return A [SIGNAL-class] object.
#' @examples
#' data(SIGNAL_example)
#' obj <- createSIGNALobject(SIGNAL_example$counts, SIGNAL_example$location,
#'                           SIGNAL_example$celltype)
#' obj
#' @export
createSIGNALobject <- function(counts, location, celltype, project = "SIGNAL") {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (length(dim(counts)) != 2) stop("`counts` must be a genes x cells matrix.")
  n <- ncol(counts)
  if (is.null(rownames(counts))) rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell", seq_len(n))

  location <- as.matrix(location)
  if (nrow(location) != n || ncol(location) < 2) {
    stop("`location` must have one row per cell (column of `counts`) and 2 columns.")
  }
  location <- location[, 1:2, drop = FALSE]
  storage.mode(location) <- "double"
  dimnames(location) <- list(colnames(counts), c("x", "y"))

  X <- celltype_matrix(celltype)
  if (nrow(X) != n) stop("`celltype` must have one entry per cell.")
  rownames(X) <- colnames(counts)
  if (anyNA(location) || anyNA(X)) stop("Missing values in `location` or `celltype`.")

  methods::new("SIGNAL", project = project, counts = counts, location = location,
               celltype = X, params = list())
}

#' Filter and transform expression
#'
#' Removes lowly expressed genes (and, with `filter = "celina"`, cells with
#' low total counts), then converts each remaining gene to normal scores.
#'
#' @details
#' With `transform = "impute"` (default, as in the code used for the
#' simulation studies) each gene's counts are converted to normal scores
#' `qnorm((rank(y) - 0.5) / n)` with tied zeros broken at random, and the
#' scores are averaged over `n_imp` draws ([rank_int_impute()]); set `seed` to
#' make this reproducible. This keeps the spatial test calibrated when many
#' counts are zero. With `transform = "rank"` (as in the original SIGNAL code
#' used for the real-data analyses) tied counts, in particular the zeros, share
#' their average rank, so the result is deterministic, but the spatial test can
#' be anti-conservative. Genes whose transformed values are constant are
#' removed.
#'
#' @param object A [SIGNAL-class] object.
#' @param filter Filtering rule:
#'   * `"min_cells"` (default): remove genes with non-zero counts in fewer
#'     than `min_cells` cells. No cells are removed.
#'   * `"celina"`: CELINA's `FilteringGenesCells()` rule. Keep genes with
#'     non-zero counts in at least `floor(percentage_cell * n_cells)` cells,
#'     then keep cells whose total count over the kept genes is greater than
#'     `min_total_counts`. Removed cells are dropped from the counts,
#'     locations and cell types.
#' @param min_cells Used by `filter = "min_cells"`.
#' @param percentage_cell,min_total_counts Used by `filter = "celina"`.
#' @param transform `"impute"` (default) or `"rank"`, see Details.
#' @param n_imp Number of zero imputations averaged when `transform = "impute"`.
#' @param seed Optional random seed, used when `transform = "impute"`.
#' @param verbose Print progress messages.
#' @return The object with `gene_expression` filled in (and, for
#'   `filter = "celina"`, `counts`, `location` and `celltype` restricted to the
#'   kept cells).
#' @export
SIGNAL_preprocess <- function(object, filter = c("min_cells", "celina"),
                              min_cells = 10, percentage_cell = 0.1,
                              min_total_counts = 10, transform = c("impute", "rank"),
                              n_imp = 10, seed = NULL, verbose = TRUE) {
  filter <- match.arg(filter)
  transform <- match.arg(transform)
  counts <- object@counts
  n_expressed <- if (inherits(counts, "dgCMatrix")) {
    tabulate(counts@i[counts@x > 0] + 1L, nbins = nrow(counts))
  } else {
    rowSums(as.matrix(counts) > 0)
  }
  min_expr <- if (filter == "celina") floor(percentage_cell * ncol(counts)) else min_cells
  keep <- n_expressed >= min_expr
  if (!any(keep)) stop("No gene is expressed in at least ", min_expr, " cells.")
  if (verbose && any(!keep)) {
    message("Removed ", sum(!keep), " genes expressed in fewer than ", min_expr, " cells.")
  }
  counts <- as.matrix(counts[keep, , drop = FALSE])
  storage.mode(counts) <- "double"
  if (anyNA(counts)) stop("Missing values in `counts`.")

  if (filter == "celina") {
    keep_cells <- colSums(counts) > min_total_counts
    if (sum(keep_cells) < 3) stop("Fewer than 3 cells pass `min_total_counts`.")
    if (any(!keep_cells)) {
      if (verbose) {
        message("Removed ", sum(!keep_cells), " cells with total count <= ", min_total_counts, ".")
      }
      counts <- counts[, keep_cells, drop = FALSE]
      object@counts <- object@counts[, keep_cells, drop = FALSE]
      object@location <- object@location[keep_cells, , drop = FALSE]
      object@celltype <- object@celltype[keep_cells, , drop = FALSE]
      if (length(object@spatial_kernels) > 0) {
        object@spatial_kernels <- list()
        object@celltype_kernel <- matrix(numeric(0), 0, 0)
        warning("Cells were removed, so the kernels were cleared; run SIGNAL_kernel() again.",
                call. = FALSE)
      }
    }
  }

  if (verbose) message("Transforming ", nrow(counts), " genes (", transform, ")")
  if (transform == "rank") {
    Y <- t(apply(counts, 1, function(y) stats::qnorm((rank(y) - 0.5) / length(y))))
  } else {
    if (!is.null(seed)) set.seed(seed)
    Y <- t(apply(counts, 1, rank_int_impute, n_imp = n_imp))
  }
  dimnames(Y) <- dimnames(counts)
  constant <- apply(Y, 1, function(y) all(y == y[1]))
  if (any(constant)) {
    if (all(constant)) stop("Every gene is constant after transformation.")
    if (verbose) message("Removed ", sum(constant), " genes with constant expression.")
    Y <- Y[!constant, , drop = FALSE]
  }

  object@gene_expression <- Y
  object@params$preprocess <- list(filter = filter, min_cells = min_cells,
                                   percentage_cell = percentage_cell,
                                   min_total_counts = min_total_counts,
                                   transform = transform, n_genes = nrow(Y),
                                   n_cells = ncol(Y), n_imp = n_imp, seed = seed)
  object
}

#' Build the cell-type and spatial kernels
#'
#' @param object A [SIGNAL-class] object.
#' @param kernels Optional list of n x n spatial kernel matrices (or a single
#'   matrix). By default the spatial kernels are built from the cell
#'   coordinates, as set by `approximation`.
#' @param approximation Which spatial kernels to build when `kernels` is
#'   `NULL`:
#'   * `"auto"` (default): the 11 kernels of [build_spatial_kernels()] for up
#'     to `max_cells` cells, and one Gaussian kernel
#'     ([build_gaussian_kernel()]) for more, as CELINA recommends for large
#'     data.
#'   * `FALSE`: always the 11 kernels.
#'   * `TRUE`: always the single Gaussian kernel.
#' @param max_cells Largest number of cells for which `approximation = "auto"`
#'   uses the 11 kernels.
#' @param bandwidth,tol Bandwidth and sparsity cut-off of the single Gaussian
#'   kernel, see [build_gaussian_kernel()]. The default bandwidth is chosen
#'   from the counts with [select_bandwidth()].
#' @param celltype_distance Distance used by the cell-type kernel,
#'   `"auto"` (default; Hamming for 0/1 single-cell labels, Euclidean for
#'   proportions), `"hamming"` or `"euclidean"`; see [build_celltype_kernel()].
#' @param basis_number,standardize Passed to [build_spatial_kernels()]
#'   (`standardize` also to [build_gaussian_kernel()]).
#' @param verbose Print progress messages.
#' @return The object with `celltype_kernel` and `spatial_kernels` filled in.
#' @export
SIGNAL_kernel <- function(object, kernels = NULL, approximation = "auto", max_cells = 6000,
                          bandwidth = NULL, tol = 1e-4,
                          celltype_distance = c("auto", "hamming", "euclidean"), basis_number = 4,
                          standardize = TRUE, verbose = TRUE) {
  celltype_distance <- match.arg(celltype_distance)
  if (celltype_distance == "auto") {
    celltype_distance <- if (all(object@celltype %in% c(0, 1))) "hamming" else "euclidean"
  }
  n <- nrow(object@location)
  user_kernels <- !is.null(kernels)
  if (!user_kernels) {
    if (identical(approximation, "auto")) approximation <- n > max_cells
    if (!is.logical(approximation) || length(approximation) != 1 || is.na(approximation)) {
      stop("`approximation` must be \"auto\", TRUE or FALSE.")
    }
    if (approximation) {
      if (verbose) message("Building one Gaussian spatial kernel for ", n, " cells")
      genes <- if (length(object@gene_expression) > 0) rownames(object@gene_expression) else NULL
      K <- build_gaussian_kernel(object@location, bandwidth = bandwidth, counts = object@counts,
                                 genes = genes, tol = tol, standardize = standardize)
      bandwidth <- attr(K, "bandwidth")
      if (verbose) message("Bandwidth: ", signif(bandwidth, 4))
      attr(K, "bandwidth") <- NULL
      kernels <- list(gaussian = K)
    } else {
      if (verbose) message("Building spatial kernels for ", n, " cells")
      kernels <- build_spatial_kernels(object@location, basis_number = basis_number,
                                       standardize = standardize)
    }
  } else {
    if (!is.list(kernels)) kernels <- list(kernels)
    kernels <- lapply(kernels, as.matrix)
    if (!all(vapply(kernels, function(K) all(dim(K) == c(n, n)), logical(1)))) {
      stop("Every kernel in `kernels` must be an n x n matrix (n = number of cells).")
    }
    if (is.null(names(kernels))) names(kernels) <- paste0("kernel_", seq_along(kernels))
  }
  object@celltype_kernel <- build_celltype_kernel(object@celltype, distance = celltype_distance)
  object@spatial_kernels <- kernels
  object@params$kernel <- list(n_kernels = length(kernels), basis_number = basis_number,
                               standardize = standardize, user_kernels = user_kernels,
                               approximation = if (user_kernels) NA else approximation,
                               bandwidth = if (!user_kernels && approximation) bandwidth else NULL,
                               tol = tol, celltype_distance = celltype_distance)
  object
}

#' Test and classify spatially variable genes
#'
#' For every gene, fits the variance component model
#' \deqn{y \sim N(0, \tau_1 K_{celltype} + \tau_2 K_{spatial} +
#'   \tau_3 K_{interaction} + \sigma^2 I)}
#' once for each spatial kernel. It runs an overall test of
#' \eqn{\tau_1 = \tau_2 = \tau_3 = 0} and tests each effect on its own:
#' \describe{
#'   \item{celltype}{expression differs between cell types;}
#'   \item{spatial}{expression varies over space, shared by all cell types;}
#'   \item{ctSVG}{the spatial pattern differs between cell types.}
#' }
#' P-values are combined across spatial kernels with the Cauchy combination
#' test ([combine_pvalues()]), and genes are classified with
#' [SIGNAL_classify()].
#'
#' @details
#' The interaction kernel is the element-wise product of the mean-centred
#' cell-type and spatial kernels. All three effect kernels are shifted to zero
#' grand mean and an identity kernel models the residual. Variance components
#' are estimated by MINQUE(0). Tail probabilities use Liu's moment-matching
#' approximation or Davies' algorithm (see `method`).
#'
#' @param object A [SIGNAL-class] object after [SIGNAL_preprocess()] and
#'   [SIGNAL_kernel()].
#' @param alpha,adjust Passed to [SIGNAL_classify()].
#' @param method How tail probabilities are computed.
#'   * `"davies"` (default): Davies' algorithm on the per-gene eigenvalues for
#'     all four tests, as in the code used for the simulation studies.
#'     Accurate down to about 1e-6, below which it reports 0 (shown as about
#'     1e-16 after combining kernels). Slowest: an n x n eigendecomposition per
#'     gene, test and kernel.
#'   * `"saddlepoint"`: saddlepoint approximation (Kuonen, 1999) on the same
#'     eigenvalues, for all four tests. As fast as `"davies"` and accurate far
#'     into the tail, so very small p-values stay distinct instead of being
#'     reported as 0. Davies' value is used near the centre of the
#'     distribution, where the saddlepoint formula is unstable.
#'   * `"liu"`: Liu's approximation for all four tests, from trace moments
#'     computed once per kernel. Fastest; its cost hardly depends on the
#'     number of genes. Can be conservative for the component tests when cell
#'     types and spatial domains overlap strongly.
#'   * `"davies_liu"`: as in the original SIGNAL code used for the real-data
#'     analyses. Davies for the overall test and Liu for the cell-type,
#'     spatial and ctSVG tests.
#' @param ncores Number of cores (forked workers; on Windows only 1 core is
#'   used). How they are used is set by `parallel`.
#' @param parallel How the cores are used:
#'   * `"genes"` (default): kernels are processed one at a time and the genes
#'     are split over the workers. The kernel set-up is computed once and
#'     shared, so each worker only adds about 3 n x n matrices, and all
#'     workers stay busy for any number of kernels.
#'   * `"kernels"`: spatial kernels are processed in parallel. Each worker
#'     holds about 25 n x n matrices, so memory grows with `ncores`, and at
#'     most one core per kernel is used. Can be faster with `"liu"` or
#'     `"davies_liu"` and few genes, whose per-kernel set-up is costly.
#'
#'   With an optimised multithreaded BLAS, set its number of threads so that
#'   threads x `ncores` does not exceed the number of cores.
#' @param verbose Print progress messages.
#' @return The object with `result` and `pvalues_by_kernel` filled in.
#'   `result` has one row per gene with columns `gene`, `p_overall`,
#'   `p_celltype`, `p_spatial`, `p_ctSVG` and the classification columns
#'   described in [SIGNAL_classify()].
#' @export
SIGNAL_test <- function(object, alpha = 0.05, adjust = "BH",
                        method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = 1,
                        parallel = c("genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .check_ready(object)
  Y <- object@gene_expression
  kernels <- object@spatial_kernels
  K_celltype <- object@celltype_kernel
  cores <- .split_cores(ncores, parallel)

  per_kernel <- .map_kernels(kernels, cores$kernels, verbose, function(Ks) {
    .signal_engine(Y, .signal_kernels(K_celltype, Ks), 1:3, method, ncores = cores$genes)
  })

  tests <- c("overall", "celltype", "spatial", "ctSVG")
  pvalues_by_kernel <- lapply(seq_along(tests), function(t) {
    P <- vapply(per_kernel, function(x) x[, t], numeric(nrow(Y)))
    matrix(P, nrow = nrow(Y), dimnames = list(rownames(Y), names(kernels)))
  })
  names(pvalues_by_kernel) <- tests

  result <- data.frame(
    gene = rownames(Y),
    p_overall = combine_pvalues(pvalues_by_kernel$overall),
    p_celltype = combine_pvalues(pvalues_by_kernel$celltype),
    p_spatial = combine_pvalues(pvalues_by_kernel$spatial),
    p_ctSVG = combine_pvalues(pvalues_by_kernel$ctSVG),
    row.names = rownames(Y),
    stringsAsFactors = FALSE
  )
  object@result <- result
  object@pvalues_by_kernel <- pvalues_by_kernel
  object@params$test <- list(method = method)
  SIGNAL_classify(object, alpha = alpha, adjust = adjust)
}

#' Classify genes from SIGNAL p-values
#'
#' A gene is significant when its multiplicity-adjusted overall p-value is
#' below `alpha`. Each significant gene is then labelled with every effect
#' whose p-value is below `alpha`: `CellType`, `Spatial` and/or `ctSVG`.
#' Use this to change the decision rule without refitting.
#'
#' @param x A [SIGNAL-class] object after [SIGNAL_test()], or a data frame
#'   with columns `p_overall`, `p_celltype`, `p_spatial` and `p_ctSVG`.
#' @param alpha Significance level.
#' @param adjust Multiple testing correction for the overall test, any method
#'   accepted by [stats::p.adjust()].
#' @return The object (or data frame) with columns `padj_overall`,
#'   `significant`, `celltype_effect`, `spatial_effect`, `ctSVG_effect` and
#'   `category` added to the result. `category` is `"Not significant"`, a
#'   combination such as `"Spatial+ctSVG"`, or `"Unclassified"` when the
#'   overall test is significant but no single effect is.
#' @export
SIGNAL_classify <- function(x, alpha = 0.05, adjust = "BH") {
  is_obj <- methods::is(x, "SIGNAL")
  res <- if (is_obj) x@result else x
  if (nrow(res) == 0) stop("Run SIGNAL_test() first.")
  res$padj_overall <- stats::p.adjust(res$p_overall, method = adjust)
  sig <- !is.na(res$padj_overall) & res$padj_overall < alpha
  below <- function(p) sig & !is.na(p) & p < alpha
  res$significant <- sig
  res$celltype_effect <- below(res$p_celltype)
  res$spatial_effect <- below(res$p_spatial)
  res$ctSVG_effect <- below(res$p_ctSVG)

  flags <- cbind(CellType = res$celltype_effect, Spatial = res$spatial_effect,
                 ctSVG = res$ctSVG_effect)
  label <- apply(flags, 1, function(f) paste(colnames(flags)[f], collapse = "+"))
  label[sig & label == ""] <- "Unclassified"
  label[!sig] <- "Not significant"
  res$category <- label

  if (!is_obj) return(res)
  x@result <- res
  x@params$alpha <- alpha
  x@params$adjust <- adjust
  x
}

#' Find the cell types behind each ctSVG
#'
#' Second-stage test for genes with a cell-type-specific spatial pattern. For
#' each cell type `c`, the interaction kernel is built from that cell type's
#' indicator kernel \eqn{x_c x_c^T / C} instead of the full cell-type kernel,
#' and only the interaction component is tested. A small p-value means the
#' gene has a spatial pattern specific to cell type `c`.
#'
#' @param object A [SIGNAL-class] object after [SIGNAL_kernel()] (and usually
#'   [SIGNAL_test()]).
#' @param genes Genes to test. Defaults to the genes with a significant ctSVG
#'   effect from [SIGNAL_test()], or all genes if the test has not been run.
#' @inheritParams SIGNAL_test
#' @return The object with `celltype_result`, a genes x cell-types matrix of
#'   Cauchy-combined p-values, filled in.
#' @export
SIGNAL_celltype <- function(object, genes = NULL, method = c("davies", "saddlepoint", "liu", "davies_liu"),
                            ncores = 1, parallel = c("genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .check_ready(object)
  Y <- object@gene_expression
  cores <- .split_cores(ncores, parallel)
  if (is.null(genes)) {
    genes <- if (nrow(object@result) > 0) {
      object@result$gene[object@result$ctSVG_effect]
    } else {
      rownames(Y)
    }
  }
  missing_genes <- setdiff(genes, rownames(Y))
  if (length(missing_genes) > 0) {
    stop("Genes not found after preprocessing: ",
         paste(utils::head(missing_genes), collapse = ", "))
  }
  X <- object@celltype
  out <- matrix(NA_real_, length(genes), ncol(X), dimnames = list(genes, colnames(X)))
  if (length(genes) == 0) {
    if (verbose) message("No genes to test.")
  } else {
    Yg <- Y[genes, , drop = FALSE]
    for (c in seq_len(ncol(X))) {
      if (verbose) message("Cell type ", colnames(X)[c])
      K_c <- tcrossprod(X[, c]) / ncol(X)
      per_kernel <- .map_kernels(object@spatial_kernels, cores$kernels, FALSE, function(Ks) {
        .signal_engine(Yg, .signal_kernels(object@celltype_kernel, Ks, K_c), 3, method,
                       ncores = cores$genes)[, 2]
      })
      out[, c] <- combine_pvalues(do.call(cbind, per_kernel))
    }
  }
  object@celltype_result <- out
  object@params$celltype <- list(method = method)
  object
}

#' Run the full SIGNAL pipeline
#'
#' Convenience wrapper for [createSIGNALobject()], [SIGNAL_preprocess()],
#' [SIGNAL_kernel()] and [SIGNAL_test()].
#'
#' @inheritParams createSIGNALobject
#' @inheritParams SIGNAL_preprocess
#' @inheritParams SIGNAL_kernel
#' @inheritParams SIGNAL_test
#' @return A [SIGNAL-class] object.
#' @examples
#' data(SIGNAL_example)
#' obj <- SIGNAL(SIGNAL_example$counts, SIGNAL_example$location,
#'               SIGNAL_example$celltype, seed = 1)
#' obj
#' summary(obj)
#' head(obj@result)
#' @export
SIGNAL <- function(counts, location, celltype, kernels = NULL,
                   approximation = "auto", max_cells = 6000,
                   celltype_distance = c("auto", "hamming", "euclidean"),
                   filter = c("min_cells", "celina"), min_cells = 10,
                   percentage_cell = 0.1, min_total_counts = 10,
                   transform = c("impute", "rank"), n_imp = 10,
                   alpha = 0.05, adjust = "BH",
                   method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = 1,
                   parallel = c("genes", "kernels"), seed = NULL,
                   project = "SIGNAL", verbose = TRUE) {
  obj <- createSIGNALobject(counts, location, celltype, project = project)
  obj <- SIGNAL_preprocess(obj, filter = filter, min_cells = min_cells,
                           percentage_cell = percentage_cell,
                           min_total_counts = min_total_counts,
                           transform = transform, n_imp = n_imp,
                           seed = seed, verbose = verbose)
  obj <- SIGNAL_kernel(obj, kernels = kernels, approximation = approximation,
                       celltype_distance = celltype_distance,
                       max_cells = max_cells, verbose = verbose)
  SIGNAL_test(obj, alpha = alpha, adjust = adjust, method = method,
              ncores = ncores, parallel = parallel, verbose = verbose)
}

.check_ready <- function(object) {
  if (length(object@gene_expression) == 0) stop("Run SIGNAL_preprocess() first.")
  if (length(object@spatial_kernels) == 0) stop("Run SIGNAL_kernel() first.")
}

# Cores for the kernel level and the gene level, from `ncores` and `parallel`.
.split_cores <- function(ncores, parallel = c("genes", "kernels")) {
  parallel <- match.arg(parallel)
  if (parallel == "kernels") list(kernels = ncores, genes = 1) else list(kernels = 1, genes = ncores)
}

.map_kernels <- function(kernels, ncores, verbose, FUN) {
  run <- function(i) {
    if (verbose) message("Spatial kernel ", i, "/", length(kernels), ": ", names(kernels)[i])
    FUN(kernels[[i]])
  }
  if (ncores > 1 && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(seq_along(kernels), run, mc.cores = ncores,
                              mc.preschedule = FALSE)
    failed <- vapply(res, inherits, logical(1), what = "try-error")
    if (any(failed)) stop("Kernel computation failed: ", res[[which(failed)[1]]])
    res
  } else {
    lapply(seq_along(kernels), run)
  }
}
