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
#' Hamming kernel between the cell-type profiles of every pair of cells:
#' \eqn{K_{ij} = \exp(-\gamma d_H(x_i, x_j))}, where \eqn{d_H} counts the
#' cell-type columns on which cells \eqn{i} and \eqn{j} differ. For one-hot
#' labels this is 1 for cells of the same type and \eqn{e^{-2\gamma}} otherwise.
#'
#' @param celltype Cell-type labels or a cells x cell-types matrix, see
#'   [celltype_matrix()].
#' @param gamma Kernel scale.
#' @return An n x n kernel matrix.
#' @examples
#' build_celltype_kernel(c("A", "B", "A"))
#' @export
build_celltype_kernel <- function(celltype, gamma = 1) {
  X <- celltype_matrix(celltype)
  n <- nrow(X)
  dH <- matrix(0, n, n)
  for (k in seq_len(ncol(X))) {
    dH <- dH + outer(X[, k], X[, k], "!=")
  }
  exp(-gamma * dH)
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
.signal_kernels <- function(K_celltype, K_spatial, K_interaction_celltype = K_celltype) {
  K_int <- (K_interaction_celltype - mean(K_interaction_celltype)) *
    (K_spatial - mean(K_spatial))
  list(
    K_celltype - mean(K_celltype),
    K_spatial - mean(K_spatial),
    K_int - mean(K_int),
    diag(nrow(K_spatial))
  )
}

utils::globalVariables(c("x", "y"))
