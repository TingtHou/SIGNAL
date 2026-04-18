#####################################################################
# Package: SIGNAL
# Version: 1.0.0
# Date : 2026-04-17
######################################################################

##########################################################
#       SIGNAL construct object functions and preprocessing related functions   #
##########################################################
#' @title Each SIGNAL object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot celltype_mat The raw celltype indicator or proportion matrix.
#' Rows are celltypes, columns are spots/cells.
#' @slot gene_expression_mat Original expression matrix. By default we use
#' normalizeCounts normalization in Scater R package.
#' @slot genes_list Different tested gene lists for cell types.
#' @slot project Name of the project (for record keeping).
#' @slot covariates The covariates in experiments (if any covariates is included).
#' @slot location Cell/spot spatial coordinates to compute the kernel matrix.
#' @slot kernelmat The kernel matrix for spatial relationship between locations.
#' @slot kernelmat_approx_U The calculated approximated kernel matrix for spatial relationship between locations.
#' @slot bandwidth The bandwidth in Gaussian kernel, users can also specify their preferred bandwidth.
#' @slot result The summary of cell type specific SE genes results for different cell types.
#' @slot approximation logical value to indicate which type of algorithms to be used in the analysis
#' @export

setClass("SIGNAL", slots = list(
  cell_type_compositions = "ANY",
  gene_expression = "ANY",
  genes_list = "ANY",
  project = "character",
  covariates = "ANY",
  location = "matrix",
  kernelmat = "ANY",
  kernelmat_approx_U = "ANY",
  bandwidth = "numeric",
  cell_type_compositions_kernel = "ANY",
  result = "ANY",
  original_location="ANY",
  original_gene_expression="ANY",
  original_cell_type_compositions="ANY",
  approximation = "logical"
))

#' @title Create the SIGNAL object.
#' @param cell_type_compositions The raw celltype indicator or proportion matrix. Rows are celltypes, columns are spots/cells.
#' @param gene_expression Original expression matrix. By default we use
#' normalizeCounts normalization in Scater R package.
#' @param location Spatial location matrix (matrix), the dimension is n x d, n is the number of locations, d is dimensin of spatial coordinates,
#' e.g. d = 2 for locations on 2D space. The rownames of locations and the colnames of count matrix should be matched.
#' @param covariates The covariates in experiments (matrix, if any covariates included), n x q,
#' n is the number of locations, q is the number of covariates.
#' The rownames of covariates and the rownames of locations should be matched.
#' @param project Name of the project (User defined).
#' @return Returns Celina object.
#'
#' @export
Create_SIGNAL_Object <- function(cell_type_compositions = NULL,
                                 gene_expression = NULL,
                                 location = NULL,
                                 covariates = NULL,
                                 project = "SIGNAL"){

  ## check dimension
  if (nrow(cell_type_compositions) != nrow(location) | ncol(gene_expression) != nrow(location)) {
    stop ("The number of columns in two matrices, and the number of locations should be consistent.")
  } # end fi

  ## check names
  if (any(rownames(cell_type_compositions) != rownames(location)) | any(colnames(gene_expression) != rownames(location))) {
    stop ("The cells/spots names should match between celltype matrix, gene expression matrix and location matrix")
  }# end fi

  ## inheriting
  object <- methods::new(
    Class = "SIGNAL",
    cell_type_compositions = cell_type_compositions,
    gene_expression = gene_expression,
    location = location,
    project = project
  )

  if (!is.null(covariates)) {
    ## check data order should consistent
    if (!identical(rownames(covariates), rownames(location))) {
      stop("The row names of covariates and row names of location should be should be matched (covariates -- n locations x k covariates; location -- n locations x d dimension).")
    }
    q <- dim(covariates)[2]
    n_covariates <- dim(covariates)[1]
    ## remove the intercept if added by user, later intercept will add automatically
    if (length(unique(covariates[, 1])) == 1){
      covariates <- covariates[, -1]
      q <- q - 1
    }
    object@covariates <- as.matrix(covariates, n_covariates, q)
  } else {
    object@covariates <- NULL
  }

  object@result <- list()
  # if (scaling_method == "separate") {
  #   object@location <- scale(object@location)
  # } else if (scaling_method == "joint") {
  #   scale_vals <- max(location) - min(location)
  #   object@location[, 1] <- (location[, 1] - min(location))/scale_vals
  #   object@location[, 2] <- (location[, 2] - min(location))/scale_vals
  # }


  ## Return created object
  return(object)
}
