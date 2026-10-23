#' The SIGNAL class
#'
#' Holds the input data, intermediate quantities and results of a SIGNAL
#' analysis. Create one with [createSIGNALobject()].
#'
#' @slot project Project name.
#' @slot counts Genes x cells count matrix as supplied (base or sparse matrix).
#' @slot location Cells x 2 matrix of spatial coordinates.
#' @slot celltype Cells x cell-types indicator (or proportion) matrix.
#' @slot gene_expression Genes x cells matrix of transformed expression, set by
#'   [SIGNAL_preprocess()].
#' @slot celltype_kernel Cell-type kernel, set by [SIGNAL_kernel()].
#' @slot spatial_kernels List of spatial kernels, set by [SIGNAL_kernel()].
#' @slot result Genes x p-values data frame (`p_overall`, `p_celltype`,
#'   `p_spatial`, `p_ctSVG`), set by [SIGNAL_test()].
#' @slot pvalues_by_kernel List of genes x kernels p-value matrices for the
#'   overall, cell-type, spatial and ctSVG tests, set by [SIGNAL_test()].
#' @slot celltype_result Genes x cell-types matrix of ctSVG p-values, set by
#'   [SIGNAL_celltype()].
#' @slot params Settings used in each step.
#' @exportClass SIGNAL
setClass("SIGNAL", slots = c(
  project = "character",
  counts = "ANY",
  location = "matrix",
  celltype = "matrix",
  gene_expression = "matrix",
  celltype_kernel = "matrix",
  spatial_kernels = "list",
  result = "data.frame",
  pvalues_by_kernel = "list",
  celltype_result = "matrix",
  params = "list"
))

#' @rdname SIGNAL-class
#' @param object A SIGNAL object.
#' @exportMethod show
setMethod("show", "SIGNAL", function(object) {
  cat("SIGNAL object:", object@project, "\n")
  cat(sprintf("  %d genes x %d cells, %d cell types\n", nrow(object@counts),
              ncol(object@counts), ncol(object@celltype)))
  steps <- c(
    preprocess = length(object@gene_expression) > 0,
    kernel = length(object@spatial_kernels) > 0,
    test = nrow(object@result) > 0,
    celltype = length(object@celltype_result) > 0
  )
  cat("  Steps run:", if (any(steps)) paste(names(steps)[steps], collapse = ", ") else "none", "\n")
  if (steps[["kernel"]]) {
    cat(sprintf("  %d spatial kernels\n", length(object@spatial_kernels)))
  }
  if (steps[["test"]]) {
    cat(sprintf("  %d genes tested; p-values in @result\n", nrow(object@result)))
  }
  invisible(object)
})
