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

#' Test spatially variable genes
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
#' test ([combine_pvalues()]).
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
#' @param ncores Number of cores. Defaults to all available cores: the
#'   Slurm allocation (`SLURM_CPUS_PER_TASK`) when set, otherwise
#'   [parallel::detectCores()] (2 under `R CMD check`). The workers are
#'   forked processes that share memory ([pbmcapply::pbmclapply()]). Windows
#'   cannot fork, so there the computation runs on one core. How the cores are
#'   used is set by `parallel`. The number of BLAS threads is left to the
#'   system: with a multithreaded BLAS (e.g. OpenBLAS), set
#'   `OPENBLAS_NUM_THREADS=1` before starting R so that the workers do not
#'   compete for cores.
#' @param parallel How the cores are used:
#'   * `"genes"` (default): kernels are processed one at a time and the genes
#'     are split over the workers. The kernel set-up is computed once and
#'     shared, so each worker only adds about 3 n x n matrices, and all
#'     workers stay busy for any number of kernels.
#'   * `"kernels"`: spatial kernels are processed in parallel. Each worker
#'     holds about 25 n x n matrices, so memory grows with `ncores`, and at
#'     most one core per kernel is used. Can be faster with `"liu"` or
#'     `"davies_liu"` and few genes, whose per-kernel set-up is costly.
#' @param verbose Print progress messages and, with `ncores > 1`, a progress
#'   bar for each spatial kernel ([pbmcapply::pbmclapply()]).
#' @return The object with `result` and `pvalues_by_kernel` filled in.
#'   `result` is a data frame with one row per gene (row names) and the
#'   p-values `p_overall`, `p_celltype`, `p_spatial` and `p_ctSVG`, combined
#'   over spatial kernels and not adjusted for multiple testing.
#' @export
SIGNAL_test <- function(object,
                        method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                        parallel = c("genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .check_ready(object)
  Y <- object@gene_expression
  kernels <- object@spatial_kernels
  K_celltype <- object@celltype_kernel
  if (is.null(ncores)) ncores <- .default_cores()
  cores <- .split_cores(ncores, parallel)

  job <- .kernel_job(Y, K_celltype, 1:3, method, ncores = cores$genes, progress = verbose)
  per_kernel <- .map_kernels(kernels, cores$kernels, verbose, job)

  tests <- c("overall", "celltype", "spatial", "ctSVG")
  pvalues_by_kernel <- lapply(seq_along(tests), function(t) {
    P <- vapply(per_kernel, function(x) x[, t], numeric(nrow(Y)))
    matrix(P, nrow = nrow(Y), dimnames = list(rownames(Y), names(kernels)))
  })
  names(pvalues_by_kernel) <- tests

  object@result <- data.frame(
    p_overall = combine_pvalues(pvalues_by_kernel$overall),
    p_celltype = combine_pvalues(pvalues_by_kernel$celltype),
    p_spatial = combine_pvalues(pvalues_by_kernel$spatial),
    p_ctSVG = combine_pvalues(pvalues_by_kernel$ctSVG),
    row.names = rownames(Y)
  )
  object@pvalues_by_kernel <- pvalues_by_kernel
  object@params$test <- list(method = method)
  object
}

#' Find the cell types behind each ctSVG
#'
#' Second-stage test for genes with a cell-type-specific spatial pattern. For
#' each cell type `c`, the interaction kernel is built from that cell type's
#' indicator kernel \eqn{x_c x_c^T / C} instead of the full cell-type kernel,
#' and only the interaction component is tested. A small p-value means the
#' gene has a spatial pattern specific to cell type `c`.
#'
#' @param object A [SIGNAL-class] object after [SIGNAL_kernel()].
#' @param genes Genes to test, usually those with a small ctSVG p-value from
#'   [SIGNAL_test()]. Defaults to all genes.
#' @inheritParams SIGNAL_test
#' @return The object with `celltype_result`, a genes x cell-types matrix of
#'   Cauchy-combined p-values, filled in.
#' @export
SIGNAL_celltype <- function(object, genes = NULL, method = c("davies", "saddlepoint", "liu", "davies_liu"),
                            ncores = NULL, parallel = c("genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .check_ready(object)
  Y <- object@gene_expression
  if (is.null(ncores)) ncores <- .default_cores()
  cores <- .split_cores(ncores, parallel)
  if (is.null(genes)) genes <- rownames(Y)
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
      job <- .kernel_job(Yg, object@celltype_kernel, 3, method, ncores = cores$genes,
                         K_interaction = K_c, column = 2)
      per_kernel <- .map_kernels(object@spatial_kernels, cores$kernels, FALSE, job)
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
#'               SIGNAL_example$celltype, seed = 1, ncores = 2)
#' obj
#' head(obj@result)
#' @export
SIGNAL <- function(counts, location, celltype, kernels = NULL,
                   approximation = "auto", max_cells = 6000,
                   celltype_distance = c("auto", "hamming", "euclidean"),
                   filter = c("min_cells", "celina"), min_cells = 10,
                   percentage_cell = 0.1, min_total_counts = 10,
                   transform = c("impute", "rank"), n_imp = 10,
                   method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
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
  SIGNAL_test(obj, method = method,
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

# All available cores: the Slurm allocation if set, otherwise the machine's
# cores; 2 under R CMD check (_R_CHECK_LIMIT_CORES_).
.default_cores <- function() {
  if (!tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false")) %in% c("", "false")) return(2L)
  slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (!is.na(slurm) && slurm > 0) return(slurm)
  max(1L, parallel::detectCores(), na.rm = TRUE)
}

# P-values for one spatial kernel.
.kernel_job <- function(Y, K_celltype, test_index, method, ncores = 1, progress = FALSE,
                        K_interaction = NULL, column = NULL) {
  force(Y); force(K_celltype)
  function(Ks) {
    K_int <- if (is.null(K_interaction)) K_celltype else K_interaction
    P <- .signal_engine(Y, .signal_kernels(K_celltype, Ks, K_int), test_index, method,
                        ncores = ncores, progress = progress)
    if (is.null(column)) P else P[, column]
  }
}

.map_kernels <- function(kernels, ncores, verbose, FUN) {
  if (ncores > 1) {
    if (verbose) message("Spatial kernels 1-", length(kernels), " on ", ncores, " workers")
    return(.par_lapply(kernels, FUN, ncores, progress = verbose, preschedule = FALSE))
  }
  lapply(seq_along(kernels), function(i) {
    if (verbose) message("Spatial kernel ", i, "/", length(kernels), ": ", names(kernels)[i])
    FUN(kernels[[i]])
  })
}
