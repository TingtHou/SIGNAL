#' Convert cell-type labels to an indicator matrix
#'
#' @param celltype A vector (character or factor) of cell-type labels, one per
#'   cell, or an already formed cells x cell-types matrix (indicators or
#'   proportions), which is returned unchanged.
#' @return A numeric matrix with one row per cell and one column per cell type.
#' @examples
#' celltype_matrix(c("A", "B", "A", "C"))
#' @export
celltype_matrix <- function(celltype) {
  if (is.matrix(celltype) || is.data.frame(celltype)) {
    X <- as.matrix(celltype)
    storage.mode(X) <- "double"
    if (is.null(colnames(X))) colnames(X) <- paste0("CT", seq_len(ncol(X)))
    return(X)
  }
  f <- factor(celltype)
  X <- vapply(levels(f), function(l) as.numeric(f == l), numeric(length(f)))
  X <- matrix(X, nrow = length(f), dimnames = list(names(celltype), levels(f)))
  X
}

#' Cell-type kernel
#'
#' Kernel between the cell-type profiles \eqn{x_i} (single-cell labels or
#' spot-level cell-type proportions) of every pair of cells,
#' \eqn{K_{ij} = \exp(-\gamma d(x_i, x_j))}, with
#' \describe{
#'   \item{`"hamming"`}{\eqn{d} the number of cell-type columns on which the
#'     two cells differ (Hamming kernel), for single-cell data where the
#'     cell-type matrix is 0/1;}
#'   \item{`"euclidean"`}{\eqn{d} the squared Euclidean distance (Gaussian
#'     kernel), for cell-type proportions, e.g. from deconvolution of spots.}
#' }
#' With `distance = "auto"` (default) the Hamming kernel is used when every
#' entry of the cell-type matrix is 0 or 1, and the Euclidean kernel
#' otherwise.
#'
#' @param celltype Cell-type labels or a cells x cell-types matrix, see
#'   [celltype_matrix()].
#' @param gamma Kernel scale.
#' @param distance `"auto"` (default), `"hamming"` or `"euclidean"`.
#' @return An n x n kernel matrix.
#' @examples
#' build_celltype_kernel(c("A", "B", "A"))
#' prop <- rbind(c(0.7, 0.3, 0), c(0.6, 0.2, 0.2), c(0, 0.1, 0.9))
#' build_celltype_kernel(prop)
#' @export
build_celltype_kernel <- function(celltype, gamma = 1, distance = c("auto", "hamming", "euclidean")) {
  distance <- match.arg(distance)
  X <- celltype_matrix(celltype)
  if (distance == "auto") distance <- if (all(X %in% c(0, 1))) "hamming" else "euclidean"
  if (distance == "euclidean") {
    d <- as.matrix(stats::dist(X))^2
  } else {
    n <- nrow(X)
    d <- matrix(0, n, n)
    for (k in seq_len(ncol(X))) {
      d <- d + outer(X[, k], X[, k], "!=")
    }
  }
  dimnames(d) <- NULL
  exp(-gamma * d)
}

#' Spatial kernels
#'
#' Builds the default panel of 11 spatial kernels used by SIGNAL from the
#' cell coordinates. With \eqn{d_{ij}} the Euclidean distance between cells
#' \eqn{i} and \eqn{j} and scales \eqn{l_1, \dots, l_5} set to the 20%, 40%,
#' ..., 100% quantiles of the absolute (standardized) coordinates, the panel is
#' \describe{
#'   \item{`gaussian_s`}{\eqn{\exp(-d_{ij}^2 / (2 l_s^2))}, s = 1, ..., 5;}
#'   \item{`matern_s`}{Matern kernel with smoothness 3/2,
#'     \eqn{(1 + l_s d_{ij}) \exp(-l_s d_{ij})}, s = 1, ..., 5;}
#'   \item{`spline`}{projection (hat) matrix \eqn{X (X^T X)^{-1} X^T} of a
#'     thin-plate regression spline basis \eqn{X} in the two coordinates,
#'     built with [mgcv::smoothCon()].}
#' }
#'
#' @param location An n x 2 matrix or data frame of cell coordinates.
#' @param basis_number Basis dimension of the thin-plate spline kernel.
#' @param standardize If `TRUE` (default), each coordinate is centred and
#'   scaled to unit variance before the kernels are built.
#' @return A named list of n x n kernel matrices.
#' @examples
#' loc <- matrix(runif(40), ncol = 2)
#' length(build_spatial_kernels(loc))
#' @export
build_spatial_kernels <- function(location, basis_number = 4, standardize = TRUE) {
  location <- as.matrix(location)[, 1:2, drop = FALSE]
  storage.mode(location) <- "double"
  if (standardize) location <- scale(location)
  ED <- as.matrix(stats::dist(location))
  lrang <- stats::quantile(abs(location), probs = seq(0.2, 1, by = 0.2))

  kernels <- list()
  for (i in seq_along(lrang)) {
    kernels[[paste0("gaussian_", i)]] <- exp(-ED^2 / (2 * lrang[[i]]^2))
    kernels[[paste0("matern_", i)]] <- .matern(ED, alpha = lrang[[i]], nu = 1.5)
  }

  centered <- sweep(location, 2, colMeans(location), "-")
  centered <- as.data.frame(centered / stats::sd(centered))
  colnames(centered) <- c("x", "y")
  sm <- mgcv::smoothCon(mgcv::s(x, y, k = basis_number, fx = TRUE, bs = "tp"),
                        data = centered)[[1]]
  mm <- sm$X
  kernels[["spline"]] <- mm %*% solve(crossprod(mm), t(mm))

  lapply(kernels, function(K) {
    dimnames(K) <- NULL
    K
  })
}

#' Single Gaussian spatial kernel for large data
#'
#' Builds one Gaussian kernel \eqn{\exp(-d_{ij}^2 / h)} from the standardized
#' coordinates, as in CELINA's `Calculate_Kernel(approximation = TRUE,
#' sparseKernel = TRUE)`. Entries below `tol` are set to 0. SIGNAL uses this
#' kernel instead of the panel of 11 kernels for large data sets (see
#' [SIGNAL_kernel()]), which cuts run time and memory about 11-fold.
#'
#' @param location An n x 2 matrix or data frame of cell coordinates.
#' @param bandwidth Bandwidth \eqn{h}. By default chosen from `counts` with
#'   [select_bandwidth()].
#' @param counts Genes x cells count matrix, used to choose the bandwidth when
#'   `bandwidth` is `NULL`.
#' @param genes Passed to [select_bandwidth()].
#' @param tol Entries of the kernel below `tol` are set to 0.
#' @param standardize If `TRUE` (default), each coordinate is centred and
#'   scaled to unit variance first.
#' @return An n x n kernel matrix with attribute `bandwidth`.
#' @examples
#' data(SIGNAL_example)
#' K <- build_gaussian_kernel(SIGNAL_example$location, counts = SIGNAL_example$counts)
#' attr(K, "bandwidth")
#' @export
build_gaussian_kernel <- function(location, bandwidth = NULL, counts = NULL, genes = NULL,
                                  tol = 1e-4, standardize = TRUE) {
  if (is.null(bandwidth)) {
    if (is.null(counts)) stop("Give `bandwidth` or `counts` to choose it from.")
    bandwidth <- select_bandwidth(counts, genes)
  }
  location <- as.matrix(location)[, 1:2, drop = FALSE]
  storage.mode(location) <- "double"
  if (standardize) location <- scale(location)
  K <- exp(-as.matrix(stats::dist(location))^2 / bandwidth)
  K[K < tol] <- 0
  dimnames(K) <- NULL
  attr(K, "bandwidth") <- bandwidth
  K
}

#' Kernel bandwidth from expression
#'
#' CELINA's default bandwidth (`bandwidth_select(method = "Silverman")`): the
#' median over genes of Silverman's rule of thumb ([stats::bw.nrd0()]) applied
#' to log-normalized expression, \eqn{\log_2(y / s + 1)} with library size
#' factors \eqn{s} scaled to mean 1 (as `scater::normalizeCounts()`). As in
#' CELINA, the size factors use all genes and the median is taken over the
#' tested genes.
#'
#' @param counts Genes x cells count matrix (base or sparse).
#' @param genes Genes (row names or indices of `counts`) over which the median
#'   is taken. Defaults to all genes.
#' @return The bandwidth, a single number.
#' @export
select_bandwidth <- function(counts, genes = NULL) {
  lib <- if (inherits(counts, "Matrix")) Matrix::colSums(counts) else colSums(as.matrix(counts))
  sf <- lib / mean(lib)
  if (!is.null(genes)) counts <- counts[genes, , drop = FALSE]
  counts_t <- if (inherits(counts, "Matrix")) Matrix::t(counts) else t(as.matrix(counts))
  bw <- vapply(seq_len(ncol(counts_t)), function(g) {
    y <- log2(as.numeric(counts_t[, g]) / sf + 1)
    tryCatch(stats::bw.nrd0(y), error = function(e) NA_real_)
  }, numeric(1))
  stats::median(bw, na.rm = TRUE)
}

# Matern covariance with smoothness nu on the scaled distance alpha * d;
# for nu = 1.5 this equals (1 + alpha * d) * exp(-alpha * d).
.matern <- function(d, alpha, nu) {
  d <- d * alpha
  d[d == 0] <- 1e-10
  (d^nu) * besselK(d, nu) / (2^(nu - 1) * gamma(nu))
}

# Assemble the kernel list for one spatial kernel: cell type, spatial,
# interaction and identity (residual). The three effect kernels are shifted to
# zero grand mean.
.signal_kernels <- function(K_celltype, K_spatial) {
  K_int <- (K_celltype - mean(K_celltype)) * (K_spatial - mean(K_spatial))
  list(
    K_celltype - mean(K_celltype),
    K_spatial - mean(K_spatial),
    K_int - mean(K_int),
    diag(nrow(K_spatial))
  )
}

utils::globalVariables(c("x", "y"))
