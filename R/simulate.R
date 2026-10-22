#' Simulate a small single-cell spatial data set
#'
#' Cells are placed uniformly on the unit square and assigned one of
#' `n_celltypes` cell types at random. Counts are negative binomial with mean 1
#' times a fold change that depends on the gene type:
#' \describe{
#'   \item{null}{no effect;}
#'   \item{celltype}{higher in cell type 1 everywhere;}
#'   \item{spatial}{higher in a central hotspot for all cell types;}
#'   \item{ctSVG}{higher in the hotspot for cell type 1 only.}
#' }
#'
#' @param n_cells Number of cells.
#' @param n_genes Named vector giving the number of genes of each type.
#' @param n_celltypes Number of cell types.
#' @param fold Fold change of the effect.
#' @param theta Negative binomial size (dispersion) parameter.
#' @param seed Optional random seed.
#' @return A list with `counts` (genes x cells), `location` (cells x 2),
#'   `celltype` (factor) and `gene_type` (named by gene).
#' @examples
#' sim <- simulate_signal_data(n_cells = 100, seed = 1)
#' str(sim, max.level = 1)
#' @export
simulate_signal_data <- function(n_cells = 300,
                                 n_genes = c(null = 40, celltype = 10, spatial = 10, ctSVG = 10),
                                 n_celltypes = 4, fold = 3, theta = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  location <- cbind(x = stats::runif(n_cells), y = stats::runif(n_cells))
  celltype <- factor(paste0("CT", sample(seq_len(n_celltypes), n_cells, replace = TRUE)),
                     levels = paste0("CT", seq_len(n_celltypes)))
  hotspot <- (location[, "x"] - 0.5)^2 + (location[, "y"] - 0.5)^2 < 0.3^2
  is_ct1 <- celltype == "CT1"

  effect <- list(
    null = rep(FALSE, n_cells),
    celltype = is_ct1,
    spatial = hotspot,
    ctSVG = hotspot & is_ct1
  )
  gene_type <- rep(names(n_genes), n_genes)
  names(gene_type) <- paste0("gene", seq_along(gene_type))
  counts <- t(vapply(gene_type, function(type) {
    mu <- ifelse(effect[[type]], fold, 1)
    stats::rnbinom(n_cells, size = theta, mu = mu)
  }, numeric(n_cells)))
  colnames(counts) <- rownames(location) <- paste0("cell", seq_len(n_cells))

  list(counts = counts, location = location, celltype = celltype,
       gene_type = gene_type)
}
