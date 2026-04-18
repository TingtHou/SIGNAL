##########################################################
#       SIGNAL main model related functions              #
##########################################################

#' @title Build the kernel list for one spatial kernel.
#' @description Constructs the centered kernel list \code{[K_CellType, K_Spatial,
#' K_Interaction]} for a single spatial kernel, where the interaction kernel is
#' the Hadamard product of the centered cell-type and spatial kernels. Each kernel
#' is then mean-centered over its entries.
#' @param Kspatial A single spatial kernel matrix.
#' @param KCellType The cell-type kernel matrix.
#' @return A list of three centered kernel matrices.
build_kernel_list <- function(Kspatial, KCellType) {
  Kinteraction <- (KCellType-mean(KCellType))*(Kspatial-mean(Kspatial))

  KList_In <- list(KCellType,Kspatial,Kinteraction)


  # Mean-center each kernel (no scaling)
  KList_In<-lapply(KList_In, function(K) {
    matrix(scale(as.vector(K), scale = FALSE),
           nrow = nrow(K), ncol = ncol(K))
  })
  return(KList_In)
}


#' @title Test variance components using multiple spatial kernels.
#' @description For each spatial kernel in \code{object@kernelmat} (or the single
#' approximation kernel), runs MINQUE-based variance component testing for the
#' cell-type, spatial, and interaction components. When multiple kernels are used,
#' p-values are combined across kernels via \code{CombinePValues}.
#'
#' @param object A SIGNAL object with \code{@kernelmat} (or \code{@kernelmat_approx_U}),
#'   \code{@gene_expression_mat}, \code{@count_use}, and \code{@KCellType} populated.
#' @param kernel_mat Optional list of kernel matrices to override \code{object@kernelmat}.
#' @param num_cores Number of cores for parallel computation.
#' @return The SIGNAL object with \code{@result} populated: a gene-by-3 matrix of
#'   p-values for the CellType, Spatial, and Interaction components.
#' @export
RunTesting <- function(object, kernel_mat = NULL, num_cores = 1) {

  if (!is.null(kernel_mat)) {
    object@kernelmat <- kernel_mat
  }

  KCellType <- object@cell_type_compositions_kernel
  count_use <- object@gene_expression

  component_names <- c("Overall","CellType", "Spatial", "ctSVG")

  # Pick the list of spatial kernels to iterate over
  spatial_kernels <- if (object@approximation) {
    list(object@kernelmat_approx_U)
  } else {
    object@kernelmat
  }
  per_kernel_results<-list()
  for (i in 1:length(spatial_kernels)) {
    Kspatial<-spatial_kernels[[i]]
    KList_In <- build_kernel_list(Kspatial, KCellType)

    # Component positions match KList_In order
    TestingID <- list(CellType = 1, Spatial = 2, Interaction = 3)

    pvalueKNN <- KNN_Parallel(
      scaled_counts = count_use,
      KList_In      = KList_In,
      ncore         = num_cores,
      TestingID     = TestingID
    )
    rownames(pvalueKNN)<-rownames(count_use)[pvalueKNN[,1]]
    per_kernel_results[[i]]<-pvalueKNN[,-1]
  }

  # Combine results across spatial kernels
  if (length(per_kernel_results) == 1L) {
    result_mat <- per_kernel_results[[1]]
    common_rows<-rownames(result_mat)
  } else {
    # Stack: for each component, build a gene x n_kernels matrix, then combine
    common_rows <- Reduce(intersect, lapply(per_kernel_results, rownames))
    n_components <- length(component_names)
    result_mat <- vapply(seq_len(n_components), function(c) {
      pmat <- sapply(per_kernel_results, function(mat) mat[common_rows, c])
      CombinePValues(pmat)
    }, numeric(nrow(count_use)))

  }

  rownames(result_mat) <- common_rows
  colnames(result_mat) <- component_names
  object@result <- result_mat

  return(object)
}
