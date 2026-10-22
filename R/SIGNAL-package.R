#' SIGNAL: Spatial varIable Gene detectioN And cLassification
#'
#' Detects spatially variable genes in single-cell resolution spatial
#' transcriptomics and classifies each one by whether its variation comes from
#' cell type, from space, or from a cell-type-specific spatial pattern (ctSVG).
#'
#' A typical analysis:
#' ```
#' obj <- createSIGNALobject(counts, location, celltype)
#' obj <- SIGNAL_preprocess(obj)
#' obj <- SIGNAL_kernel(obj)
#' obj <- SIGNAL_test(obj)
#' obj <- SIGNAL_celltype(obj)
#' ```
#' or [SIGNAL()] for the first four steps in one call.
#'
#' @keywords internal
#' @importFrom methods setClass setMethod new is show
"_PACKAGE"
