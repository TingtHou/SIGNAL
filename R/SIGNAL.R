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
#' With `transform = "rank"` (default) each gene's counts are converted to
#' normal scores `qnorm((rank(y) - 0.5) / n)`; tied counts share their average
#' rank, so the result is deterministic.
#'
#' **Sparse data with many zeros:** all zero counts of a gene then share one
#' score, and this large block of tied values can make the spatial test
#' anti-conservative. For such data use `transform = "impute"`: tied zeros are
#' broken at random and the scores are averaged over `n_imp` draws
#' ([rank_int_impute()]); set `seed` to make this reproducible.
#'
#' Genes whose transformed values are constant are removed.
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
#' @param transform `"rank"` (default) or `"impute"` (recommended for sparse data with
#'   many zeros), see Details.
#' @param n_imp Number of zero imputations averaged when `transform = "impute"`.
#' @param seed Optional random seed, used when `transform = "impute"`.
#' @param verbose Print progress messages.
#' @return The object with `gene_expression` filled in (and, for
#'   `filter = "celina"`, `counts`, `location` and `celltype` restricted to the
#'   kept cells).
#' @export
SIGNAL_preprocess <- function(object, filter = c("min_cells", "celina"),
                              min_cells = 10, percentage_cell = 0.1,
                              min_total_counts = 10, transform = c("rank", "impute"),
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
  if (transform == "impute" && !is.null(seed)) set.seed(seed)
  Y <- .transform_counts(counts, transform, n_imp)
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

# Normal scores of each row of a genes x cells count matrix.
.transform_counts <- function(counts, transform, n_imp) {
  Y <- if (transform == "rank") {
    t(apply(counts, 1, function(y) stats::qnorm((rank(y) - 0.5) / length(y))))
  } else {
    t(apply(counts, 1, rank_int_impute, n_imp = n_imp))
  }
  dimnames(Y) <- dimnames(counts)
  Y
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
#' @param ncores Number of BLAS threads used to build the kernels (with a
#'   multithreaded BLAS); the previous setting is restored afterwards. Defaults
#'   to all available cores, as in [SIGNAL_test()].
#' @param verbose Print progress messages.
#' @return The object with `celltype_kernel` and `spatial_kernels` filled in.
#' @export
SIGNAL_kernel <- function(object, kernels = NULL, approximation = "auto", max_cells = 6000,
                          bandwidth = NULL, tol = 1e-4,
                          celltype_distance = c("auto", "hamming", "euclidean"), basis_number = 4,
                          standardize = TRUE, ncores = NULL, verbose = TRUE) {
  celltype_distance <- match.arg(celltype_distance)
  if (is.null(ncores)) ncores <- .default_cores()
  if (ncores > 1) {
    blas <- .get_blas_threads()
    on.exit(.set_blas_threads(blas))
    .set_blas_threads(ncores)
  }
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
#'   [parallel::detectCores()] (2 under `R CMD check`). How the cores are
#'   used is set by `parallel`. With a multithreaded BLAS (e.g. OpenBLAS,
#'   MKL), the threads are set by SIGNAL (package RhpcBLASctl) and restored
#'   afterwards:
#'   * each kernel's set-up runs on `ncores` BLAS threads;
#'   * the genes are split over `min(ncores, number of genes)` forked workers
#'     that share memory ([pbmcapply::pbmclapply()]), one BLAS thread each;
#'     each worker needs about 7 n x n matrices (about 5 GB at n = 9,000), so
#'     lower `ncores` if memory is short;
#'   * a single gene, or Windows (which cannot fork), runs on `ncores` BLAS
#'     threads.
#'
#'   So `OPENBLAS_NUM_THREADS=1` (or `OMP_NUM_THREADS=1`) can be kept in the
#'   environment for code outside SIGNAL; [SIGNAL_kernel()] also builds the
#'   kernels on `ncores` BLAS threads.
#' @param parallel How the cores are used:
#'   * `"auto"` (default): whichever of the two below is expected to be
#'     faster, from the number of kernels, genes and cores and the cost of
#'     `method`'s per-kernel set-up; `"kernels"` only if its workers fit in
#'     the available memory.
#'   * `"genes"`: kernels are processed one at a time and the genes are split
#'     over the workers. The kernel set-up is computed once and shared, so
#'     each worker only adds about 3 n x n matrices, and all workers stay busy
#'     for any number of kernels.
#'   * `"kernels"`: spatial kernels are processed in parallel. Each worker
#'     does its own set-up and holds about 20 (`"davies"`) to 50 (`"liu"`,
#'     `"davies_liu"`) n x n matrices, and at most one core per kernel is
#'     used. Faster with `"liu"` or `"davies_liu"` and few genes, whose
#'     per-kernel set-up is costly.
#' @param verbose Print progress messages and, with `ncores > 1`, a progress
#'   bar for each spatial kernel ([pbmcapply::pbmclapply()]).
#' @return The object with `result` and `pvalues_by_kernel` filled in.
#'   `result` is a data frame with one row per gene (row names) and the
#'   p-values `p_overall`, `p_celltype`, `p_spatial` and `p_ctSVG`, combined
#'   over spatial kernels and not adjusted for multiple testing.
#' @seealso [SIGNAL_overall()] and [SIGNAL_component()] to run the two steps
#'   separately, e.g. the component tests only for the genes detected by the
#'   overall test.
#' @export
SIGNAL_test <- function(object,
                        method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                        parallel = c("auto", "genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .run_tests(object, genes = NULL, overall = TRUE, components = TRUE, method = method,
             ncores = ncores, parallel = parallel, verbose = verbose)
}

#' Step 1: the overall test
#'
#' Tests every gene for any cell-type or spatial variation, the overall test
#' of \eqn{\tau_1 = \tau_2 = \tau_3 = 0} in the model of [SIGNAL_test()].
#' Genes detected here can then be classified with [SIGNAL_component()].
#'
#' @inheritParams SIGNAL_test
#' @return The object with `p_overall` in `result` (one row per gene) and
#'   `pvalues_by_kernel$overall` filled in; component p-values already in
#'   `result` are kept.
#' @examples
#' data(SIGNAL_example)
#' obj <- createSIGNALobject(SIGNAL_example$counts, SIGNAL_example$location,
#'                           SIGNAL_example$celltype)
#' obj <- SIGNAL_preprocess(obj, verbose = FALSE)
#' obj <- SIGNAL_kernel(obj, verbose = FALSE)
#' obj <- SIGNAL_overall(obj, ncores = 2, verbose = FALSE)
#' detected <- rownames(obj@result)[p.adjust(obj@result$p_overall, "BH") < 0.05]
#' obj <- SIGNAL_component(obj, genes = detected, ncores = 2, verbose = FALSE)
#' head(obj@result)
#' @export
SIGNAL_overall <- function(object,
                           method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                           parallel = c("auto", "genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .run_tests(object, genes = NULL, overall = TRUE, components = FALSE, method = method,
             ncores = ncores, parallel = parallel, verbose = verbose)
}

#' Step 2: the component tests
#'
#' Tests each effect of the model of [SIGNAL_test()] on its own, for the genes
#' in `genes`: the cell-type test (\eqn{\tau_1}), the spatial test
#' (\eqn{\tau_2}) and the ctSVG test (\eqn{\tau_3}). Usually run on the genes
#' detected by the overall test ([SIGNAL_overall()]), which saves most of the
#' run time when only some genes are detected.
#'
#' @inheritParams SIGNAL_test
#' @param genes Genes to test. Defaults to all genes.
#' @return The object with `p_celltype`, `p_spatial` and `p_ctSVG` in
#'   `result` filled in for `genes` (`NA` for other genes) and the matching
#'   `pvalues_by_kernel` entries; `p_overall` is kept.
#' @export
SIGNAL_component <- function(object, genes = NULL,
                             method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                             parallel = c("auto", "genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  .run_tests(object, genes = genes, overall = FALSE, components = TRUE, method = method,
             ncores = ncores, parallel = parallel, verbose = verbose)
}

# Run the overall and/or component tests for `genes` (default all) over all
# spatial kernels, and write the Cauchy-combined p-values into the object.
# Columns of tests not run, and rows of genes not tested, keep their values
# (NA if never run).
.run_tests <- function(object, genes, overall, components, method, ncores, parallel, verbose) {
  .check_ready(object)
  Y <- object@gene_expression
  all_genes <- rownames(Y)
  if (is.null(genes)) genes <- all_genes
  missing_genes <- setdiff(genes, all_genes)
  if (length(missing_genes) > 0) {
    stop("Genes not found after preprocessing: ", paste(utils::head(missing_genes), collapse = ", "))
  }
  tests <- c("overall", "celltype", "spatial", "ctSVG")
  cols <- c("p_overall", "p_celltype", "p_spatial", "p_ctSVG")

  # Start from the existing results, or empty ones
  res <- object@result
  if (nrow(res) == 0 || !identical(rownames(res), all_genes)) {
    res <- as.data.frame(matrix(NA_real_, length(all_genes), 4, dimnames = list(all_genes, cols)))
  }
  kernels <- object@spatial_kernels
  pk <- object@pvalues_by_kernel
  if (!all(tests %in% names(pk)) || !identical(rownames(pk$overall), all_genes)) {
    pk <- stats::setNames(lapply(tests, function(t) {
      matrix(NA_real_, length(all_genes), length(kernels), dimnames = list(all_genes, names(kernels)))
    }), tests)
  }
  if (length(genes) == 0) {
    if (verbose) message("No genes to test.")
    return(object)
  }

  if (is.null(ncores)) ncores <- .default_cores()
  test_index <- if (components) 1:3 else integer(0)
  cores <- .split_cores(ncores, parallel, length(kernels), length(genes), method, ncol(Y))
  if (verbose) {
    message(if (overall && components) "Overall and component tests" else if (overall) "Overall test"
            else "Component tests", " for ", length(genes), " genes")
  }
  job <- .kernel_job(Y[genes, , drop = FALSE], object@celltype_kernel, test_index, method,
                     ncores = cores$genes, progress = verbose, overall = overall)
  per_kernel <- .map_kernels(kernels, cores$kernels, verbose, job)

  run <- c(overall = overall, celltype = components, spatial = components, ctSVG = components)
  for (t in which(run)) {
    P <- vapply(per_kernel, function(x) x[, t], numeric(length(genes)))
    P <- matrix(P, nrow = length(genes))
    pk[[tests[t]]][genes, ] <- P
    res[genes, cols[t]] <- combine_pvalues(P)
  }
  object@result <- res
  object@pvalues_by_kernel <- pk
  object@params$test <- list(method = method)
  object
}

#' Find the cell types behind each ctSVG
#'
#' Tests, for each cell type, whether a gene has a spatial pattern in that
#' cell type. Usually run on the genes with a small ctSVG p-value from
#' [SIGNAL_component()]. The test depends on the cell-type input:
#' \describe{
#'   \item{Single-cell data (0/1 cell-type labels)}{Only the cells of type
#'     `c` are used. Their expression is transformed again within the subset
#'     (with the settings of [SIGNAL_preprocess()]), the spatial kernels are
#'     built from their coordinates (with the settings of [SIGNAL_kernel()])
#'     and mean-centred, and the model
#'     \eqn{y \sim N(0, \tau K_{spatial} + \sigma^2 I)} is tested for
#'     \eqn{\tau = 0}.}
#'   \item{Spot-level data (cell-type proportions)}{All spots are used. The
#'     interaction kernel of [SIGNAL_test()] is built from cell type `c`
#'     alone, \eqn{(x_c x_c^T) \circ K_{spatial}} with both factors
#'     mean-centred, where \eqn{x_c} is the proportion of `c` in each spot.
#'     The model keeps the cell-type and spatial kernels, and the
#'     interaction component is tested.}
#' }
#' P-values are combined across spatial kernels with [combine_pvalues()].
#'
#' @param object A [SIGNAL-class] object after [SIGNAL_kernel()].
#' @param genes Genes to test, usually the ctSVGs. Defaults to all genes.
#' @param level `"auto"` (default; `"single_cell"` when every cell-type
#'   value is 0 or 1, otherwise `"spot"`), `"single_cell"` or `"spot"`.
#' @param min_cells Single-cell level: cell types with fewer cells are not
#'   tested (`NA`).
#' @inheritParams SIGNAL_test
#' @return The object with `celltype_result`, a genes x cell types matrix of
#'   p-values not adjusted for multiple testing, filled in. Genes that are
#'   constant within a cell type get `NA` for it.
#' @export
SIGNAL_celltype <- function(object, genes = NULL, level = c("auto", "single_cell", "spot"),
                            min_cells = 10,
                            method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                            parallel = c("auto", "genes", "kernels"), verbose = TRUE) {
  method <- match.arg(method)
  level <- match.arg(level)
  .check_ready(object)
  Y <- object@gene_expression
  if (is.null(genes)) genes <- rownames(Y)
  missing_genes <- setdiff(genes, rownames(Y))
  if (length(missing_genes) > 0) {
    stop("Genes not found after preprocessing: ", paste(utils::head(missing_genes), collapse = ", "))
  }
  X <- object@celltype
  if (level == "auto") level <- if (all(X %in% c(0, 1))) "single_cell" else "spot"
  if (is.null(ncores)) ncores <- .default_cores()
  out <- matrix(NA_real_, length(genes), ncol(X), dimnames = list(genes, colnames(X)))
  if (length(genes) == 0) {
    if (verbose) message("No genes to test.")
    return(object)
  }
  pp <- object@params$preprocess
  if (level == "single_cell" && identical(pp$transform, "impute") && !is.null(pp$seed)) {
    set.seed(pp$seed)
  }

  for (c in seq_len(ncol(X))) {
    if (level == "single_cell") {
      cells <- X[, c] == 1
      if (sum(cells) < min_cells) {
        if (verbose) message("Cell type ", colnames(X)[c], ": ", sum(cells), " cells, not tested")
        next
      }
      counts <- as.matrix(object@counts[genes, cells, drop = FALSE])
      storage.mode(counts) <- "double"
      Yc <- .transform_counts(counts, pp$transform, pp$n_imp)
      Yc <- Yc[apply(Yc, 1, function(y) any(y != y[1])), , drop = FALSE]
      kernels <- .with_blas_threads(ncores, .subset_spatial_kernels(object, cells))
      test <- function(Yg, Ks, n) {
        .signal_engine(Yg, list(Ks - mean(Ks), diag(nrow(Ks))), integer(0), method, ncores = n)[, 1]
      }
    } else {
      cells <- rep(TRUE, nrow(X))
      Yc <- Y[genes, , drop = FALSE]
      kernels <- object@spatial_kernels
      K_c <- tcrossprod(X[, c])
      test <- function(Yg, Ks, n) {
        .signal_engine(Yg, .signal_kernels(object@celltype_kernel, Ks, K_c), 3, method,
                       ncores = n, overall = FALSE)[, 2]
      }
    }
    if (verbose) {
      message("Cell type ", colnames(X)[c], ": ", nrow(Yc), " genes, ", sum(cells), " cells")
    }
    if (nrow(Yc) == 0) next
    cores <- .split_cores(ncores, parallel, length(kernels), nrow(Yc), method, ncol(Yc))
    per_kernel <- .map_kernels(kernels, cores$kernels, FALSE, function(Ks) test(Yc, Ks, cores$genes))
    out[rownames(Yc), c] <- combine_pvalues(matrix(unlist(per_kernel), nrow = nrow(Yc)))
  }
  object@celltype_result <- out
  object@params$celltype <- list(level = level, method = method, min_cells = min_cells)
  object
}

# Spatial kernels for a subset of cells, built as SIGNAL_kernel() built the
# full ones: from the subset's coordinates, or the rows and columns of
# user-supplied kernels. The single Gaussian kernel keeps the full data's
# bandwidth.
.subset_spatial_kernels <- function(object, cells) {
  p <- object@params$kernel
  if (isTRUE(p$user_kernels)) {
    return(lapply(object@spatial_kernels, function(K) K[cells, cells, drop = FALSE]))
  }
  location <- object@location[cells, , drop = FALSE]
  if (isTRUE(p$approximation)) {
    K <- build_gaussian_kernel(location, bandwidth = p$bandwidth, tol = p$tol,
                               standardize = p$standardize)
    attr(K, "bandwidth") <- NULL
    return(list(gaussian = K))
  }
  build_spatial_kernels(location, basis_number = p$basis_number, standardize = p$standardize)
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
                   transform = c("rank", "impute"), n_imp = 10,
                   method = c("davies", "saddlepoint", "liu", "davies_liu"), ncores = NULL,
                   parallel = c("auto", "genes", "kernels"), seed = NULL,
                   project = "SIGNAL", verbose = TRUE) {
  obj <- createSIGNALobject(counts, location, celltype, project = project)
  obj <- SIGNAL_preprocess(obj, filter = filter, min_cells = min_cells,
                           percentage_cell = percentage_cell,
                           min_total_counts = min_total_counts,
                           transform = transform, n_imp = n_imp,
                           seed = seed, verbose = verbose)
  obj <- SIGNAL_kernel(obj, kernels = kernels, approximation = approximation,
                       celltype_distance = celltype_distance,
                       max_cells = max_cells, ncores = ncores, verbose = verbose)
  SIGNAL_test(obj, method = method,
              ncores = ncores, parallel = parallel, verbose = verbose)
}

.check_ready <- function(object) {
  if (length(object@gene_expression) == 0) stop("Run SIGNAL_preprocess() first.")
  if (length(object@spatial_kernels) == 0) stop("Run SIGNAL_kernel() first.")
}

# Cores for the kernel level and the gene level, from `ncores` and `parallel`.
# "auto" compares the expected run time of the two modes, in units of one
# n x n matrix product: per-kernel set-up S and per-gene work g for `method`.
#   genes:   n_kernels * (S + n_genes * g / ncores)
#   kernels: ceiling(n_kernels / k) * (S + n_genes * g), k kernel workers
# where k = min(ncores, n_kernels) and at most the number of kernel workers
# that fit in the available memory.
.split_cores <- function(ncores, parallel = c("auto", "genes", "kernels"),
                         n_kernels = 1, n_genes = 1, method = "davies", n_cells = 1) {
  parallel <- match.arg(parallel)
  if (parallel == "genes") return(list(kernels = 1, genes = ncores))
  if (parallel == "kernels") return(list(kernels = ncores, genes = 1))
  k <- .auto_parallel(ncores, n_kernels, n_genes, method, n_cells)
  if (k > 1) list(kernels = k, genes = 1) else list(kernels = 1, genes = ncores)
}

# Number of kernel workers for "auto": > 1 for kernel-level parallelism,
# 1 for gene-level parallelism.
.auto_parallel <- function(ncores, n_kernels, n_genes, method, n_cells) {
  if (ncores <= 1 || n_kernels <= 1) return(1L)
  heavy_setup <- method %in% c("liu", "davies_liu")
  S <- if (heavy_setup) 39 else 0.5
  g <- switch(method, liu = 0.01, davies_liu = 2, 13)
  k <- min(ncores, n_kernels)
  # Each kernel worker holds about 50 (liu, davies_liu) or 20 n x n matrices
  per_worker_gb <- (if (heavy_setup) 50 else 20) * 8 * n_cells^2 / 1024^3
  avail_gb <- .available_memory_gb()
  if (!is.na(avail_gb)) k <- min(k, floor(0.9 * avail_gb / per_worker_gb))
  if (k < 2) return(1L)
  t_genes <- n_kernels * (S + n_genes * g / ncores)
  t_kernels <- ceiling(n_kernels / k) * (S + n_genes * g)
  if (t_kernels < t_genes) as.integer(k) else 1L
}

# Available memory in GB (Linux: MemAvailable in /proc/meminfo); NA if unknown.
.available_memory_gb <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  m <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  if (!length(m)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", m)) / 1024^2
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
                        overall = TRUE) {
  force(Y); force(K_celltype)
  function(Ks) {
    .signal_engine(Y, .signal_kernels(K_celltype, Ks), test_index, method,
                   ncores = ncores, progress = progress, overall = overall)
  }
}

.map_kernels <- function(kernels, ncores, verbose, FUN) {
  if (ncores > 1) {
    if (verbose) message("Spatial kernels 1-", length(kernels), " on ", ncores, " workers")
    # Kernel workers each use one BLAS thread
    blas <- .get_blas_threads()
    on.exit(.set_blas_threads(blas))
    .set_blas_threads(1)
    return(.par_lapply(kernels, FUN, ncores, progress = verbose, preschedule = FALSE))
  }
  lapply(seq_along(kernels), function(i) {
    if (verbose) message("Spatial kernel ", i, "/", length(kernels), ": ", names(kernels)[i])
    FUN(kernels[[i]])
  })
}
