#' Example single-cell spatial data set
#'
#' A simulated data set of 400 cells and 70 genes, created with
#' [simulate_signal_data()] (see `data-raw/SIGNAL_example.R`). Cells lie on
#' the unit square and belong to four cell types. A circular hotspot of radius
#' 0.3 is centred at (0.5, 0.5).
#'
#' @format A list with
#' \describe{
#'   \item{counts}{70 x 400 count matrix (genes x cells).}
#'   \item{location}{400 x 2 matrix of cell coordinates.}
#'   \item{celltype}{factor of cell types `CT1` to `CT4`.}
#'   \item{gene_type}{true type of each gene: `null` (genes 1-40), `celltype`
#'     (41-50, higher in CT1), `spatial` (51-60, higher in the hotspot) and
#'     `ctSVG` (61-70, higher in the hotspot for CT1 only). Effects are 4-fold
#'     increases in the negative binomial mean.}
#' }
"SIGNAL_example"
