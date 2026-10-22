#' Plot the spatial expression pattern of a gene
#'
#' Draws each cell at its location, coloured by the gene's relative
#' expression: counts rescaled to the range 0 to 1.
#'
#' @param object A [SIGNAL-class] object.
#' @param gene Gene name.
#' @param by_celltype If `TRUE`, draw one panel per cell type.
#' @param point_size Point size.
#' @return A `ggplot` object.
#' @examples
#' data(SIGNAL_example)
#' obj <- createSIGNALobject(SIGNAL_example$counts, SIGNAL_example$location,
#'                           SIGNAL_example$celltype)
#' SIGNAL_plotGene(obj, "gene61", by_celltype = TRUE)
#' @export
SIGNAL_plotGene <- function(object, gene, by_celltype = FALSE, point_size = 1.5) {
  if (!gene %in% rownames(object@counts)) stop("Gene ", gene, " not found.")
  y <- as.numeric(as.matrix(object@counts[gene, , drop = FALSE]))
  rng <- max(y) - min(y)
  X <- object@celltype
  df <- data.frame(
    x = object@location[, 1],
    y = object@location[, 2],
    expr = if (rng > 0) (y - min(y)) / rng else rep(0, length(y)),
    celltype = factor(colnames(X)[max.col(X, ties.method = "first")], levels = colnames(X))
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = expr)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_gradientn(
      colours = c("antiquewhite", "#35B779", "#31688E", "#440154"),
      limits = c(0, 1), name = "Relative\nexpression") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = gene) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  if (by_celltype) p <- p + ggplot2::facet_wrap(~celltype)
  p
}

utils::globalVariables(c("expr", "celltype"))
