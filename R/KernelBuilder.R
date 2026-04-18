##########################################################
#       SIGNAL calculate kernel matrices related functions             #
##########################################################
#' @title Calculate the kernel matrix.
#' @param object The SIGNAL object.
#' @param approximation Option to use the default 11 kernels combined results or 1 kernel approximated results
#' True is to proceed with approximated algorithm
#' @param bandwidth.set.by.user User could select their own bandwidth (a numeric value)
#' if the recommended bandwidth doesn't work in their dataset.
#' @param bandwidthtype The type of bandwidth to be used in Gaussian kernel,
#' "SJ" for Sheather & Jones (1991) method (usually used in small size datasets),
#' "Silverman" for Silverman's ‘rule of thumb’ method (1986) (usually used in large size datasets).
#' @param sparseKernel Select "TURE" if the user wants to use a sparse kernel matrix or "FALSE" if not.
#' It is recommended to choose sparseKernel = "TRUE" when sample size is large.
#' @param sparseKernel_tol When sparseKernel = TRUE,
#' the cut-off value when building sparse kernel matrix,
#' any element in the kernel matrix greater than sparseKernel_tol will be kept,
#' otherwise will be set to 0 to save memory.
#' @param sparseKernel_ncore When sparseKernel=TRUE,
#' the number of CPU cores to build sparse kernel matrix.
#' @return Returns SIGNAL object.
#'
#' @export
Calculate_Kernel <- function(object,
                             approximation = FALSE,
                             bandwidth.set.by.user = NULL,
                             bandwidthtype = "Silverman",
                             sparseKernel = FALSE,
                             sparseKernel_tol = 1e-10,
                             sparseKernel_ncore = 1) {

  if (approximation == F) {
    ## Default - construct 11 kernels
    kernel_res <- SIGNAL_kernel(object@location)
    object@kernelmat <- kernel_res
    object@approximation <- FALSE

    return(object)
  } else {
    ## Alternative - 1 kernels with approximation
    object@approximation <- TRUE
    sparseKernel <- TRUE
    ## Select the bandwidth for the approximation kernel
    if (is.null(bandwidth.set.by.user) & is.null(bandwidthtype)) {
      stop("Please specify a bandwidth for kernel matrix, or specify a bandwidthtype and
           provide a normalized expression (gene by location).")
    } else if(!is.null(bandwidth.set.by.user)) {
      ## User defined bandwidth
      object@bandwidth <- bandwidth.set.by.user
      cat(paste("## The bandwidth is: ", object@bandwidth, " \n"))
    } else if(!is.null(bandwidthtype)) {
      ## Silverman selected bandwidth
      object@bandwidth <- bandwidth_select(object@gene_expression_mat, method = bandwidthtype)
      cat(paste("## The bandwidth is: ", object@bandwidth, " \n"))
    }

    ## Calculate the approximated kernel matrix
    if (sparseKernel == FALSE) {
      cat(paste("## Calculating kernel matrix\n"))
      object@kernelmat <- kernel_build(kerneltype = "gaussian",
                                       location = object@location, bandwidth = object@bandwidth)

      ## SVD on the original gaussian kernel and extract the first two eigen values
      svd_results <- RSpectra::eigs_sym(object@kernelmat, k = 2)
      newkernel <- svd_results$vectors %*% diag(svd_results$values) %*% t(svd_results$vectors)
      ## Normalize the kernel
      k_norm_diag <- 1/sqrt(diag(newkernel))
      newkernel <- newkernel * (k_norm_diag %*% t(k_norm_diag))
      object@kernelmat_approx_U <- as.matrix(svd_results$vectors %*% sqrt(diag(svd_results$values)) * k_norm_diag)
    } else if (sparseKernel == TRUE) {
      cat(paste("## Calculating sparse kernel matrix\n"))
      object@kernelmat <- kernel_build_sparse(kerneltype = "gaussian",
                                              location = object@location,
                                              bandwidth = object@bandwidth,
                                              tol = sparseKernel_tol,
                                              ncores = sparseKernel_ncore)
      ## SVD on the original gaussian kernel and extract the first two eigen values
      svd_results <- RSpectra::eigs_sym(object@kernelmat, k = 2)
      newkernel <- svd_results$vectors %*% diag(svd_results$values) %*% t(svd_results$vectors)
      ## Normalize the kernel
      k_norm_diag <- 1/sqrt(diag(newkernel))
      newkernel <- newkernel * (k_norm_diag %*% t(k_norm_diag))
      object@kernelmat_approx_U <- as.matrix(svd_results$vectors %*% sqrt(diag(svd_results$values)) * k_norm_diag)
    }
    cat(paste("## Finished calculating kernel matrix.\n"))
    return(object)
  }
}



#' @title Construct the kernel matrix
#' @description This function construct the kernel matrices in a list.
#' @param location Cell/spot spatial coordinates to compute the kernel matrix.
#' @param basis_number Number of spline functions to construct non parametric covariance matrix
#' @return A list with 11 pre-calculated kernels
SIGNAL_kernel <- function(location, basis_number = 4) {
  ## Calculate the Euclidean distance matrix
  ED <- as.matrix(dist(location[ , 1:2]))

  kernels_list <- list()
  lrang <- quantile(abs(location), probs = seq(0.2, 1, by = 0.2))
  for (ikernel in c(1:5)) {
    ## Gaussian kernel
    kernel_mat <- exp(-ED ^ 2/(2 * lrang[ikernel] ^ 2))
    kernels_list <- c(kernels_list, list(kernel_mat))
    rm(kernel_mat)

    ## matern kernels
    kernel_mat <- fields::Matern(d = ED, range = 1, alpha = lrang[ikernel], smoothness = 1.5)
    kernels_list <- c(kernels_list, list(kernel_mat))
    rm(kernel_mat)
  }# end for ikernel

  ## Add in the spline basis kernels
  center_coords <- sweep(location, 2, apply(location, 2, mean), '-')
  center_coords <- as.data.frame(center_coords / sd(as.matrix(center_coords)))
  colnames(center_coords) <- c("x", "y")
  sm <- mgcv::smoothCon(mgcv::s(x, y,
                                k = basis_number, fx = T, bs = 'tp'), data = center_coords)[[1]]
  mm <- as.matrix(data.frame(sm$X))
  mm_inv <- solve(crossprod(mm, mm))
  spline_kernels <- list(mm %*% mm_inv %*% t(mm))

  ## Final output kernels
  output_kernels <- c(kernels_list, spline_kernels)
  return(output_kernels)
}

#' @title Build kernel matrix.
#' @description This function calculates kernel matrix from spatial locations.
#' @param kerneltype The type of kernel to be used, either "gaussian",
#' or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @return The kernel matrix for spatial relationship between locations.
kernel_build <- function (kerneltype = "gaussian", location, bandwidth)
{
  if (kerneltype == "gaussian") {
    K <- exp(-1 * as.matrix(stats::dist(location) ^ 2)/bandwidth)
  }
  else if (kerneltype == "cauchy") {
    K <- 1/(1 + 1 * as.matrix(stats::dist(location) ^ 2)/bandwidth)
  }
  else if (kerneltype == "quadratic") {
    ED2 <- 1 * as.matrix(stats::dist(location) ^ 2)
    K <- 1 - ED2/(ED2 + (bandwidth))
  }
  return(K)
}


#' @title Select bandwidth in Gaussian kernel.
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr A m gene by n location matrix of normalized gene expression matrix.
#' @param method The method used in bandwidth selection, "SJ" usually for small sample size data,
#' "Silverman" usually for large sample size data.
#' @return A numeric value of calculated bandwidth.
bandwidth_select <- function (expr, method = "Silverman") {
  N <- dim(expr)[2]
  if (method == "SJ") {
    bw_SJ <- c()
    for (i in 1:dim(expr)[1]) {
      tryCatch ({
        bw_SJ[i] <- stats::bw.SJ(expr[i, ], method = "dpi")
      }, error = function(e){cat("Gene", i, " :", conditionMessage(e), "\n")})
    }
    bw <- stats::median(na.omit(bw_SJ))
  } else if (method == "Silverman") {
    bw_Silverman <- c()
    for (i in 1:dim(expr)[1]) {
      tryCatch({
        bw_Silverman[i] <- stats::bw.nrd0(expr[i, ])
      }, error = function(e){cat("Gene", i, " :", conditionMessage(e), "\n")})
    }
    bw <- median(na.omit(bw_Silverman))
  }
  return(bw)
}


#' @title Build sparse kernel matrix.
#' @description This function calculates kernel matrix.
#' @param kerneltype The type of kernel to be used, either "gaussian",
#' or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @param tol A numeric value of cut-off value when building sparse kernel matrix.
#' @param ncores A integer value of number of CPU cores to use when building sparse kernel matrix.
#' @return The sparse kernel matrix for spatial relationship between locations.
kernel_build_sparse <- function(kerneltype, location, bandwidth, tol, ncores) {
  if (kerneltype == "gaussian") {
    fx_gaussian <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      line_i[-i] <- exp(-(pdist::pdist(location[i, ], location[-i, ])@dist^2)/(bandwidth))
      ind_i <- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
    }
    ## Aggregate the sparse matrix
    results <- parallel::mclapply(1:dim(location)[1], fx_gaussian, mc.cores = ncores)
    tib <- tidyr::tibble(results) %>% tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]),
                                     x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                                    dim(location)[1]))
  } else if (kerneltype == "cauchy") {
    fx_cauchy <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      line_i[-i] <- 1/(1 + (pdist::pdist(location[i, ],location[-i, ])@dist^2)/(bandwidth))
      ind_i <- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
    }
    ## Aggregate the sparse matrix
    results <- parallel::mclapply(1:dim(location)[1], fx_cauchy, mc.cores = ncores)
    tib <- tidyr::tibble(results)  %>%  tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]),
                                     x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                                    dim(location)[1]))
  } else if (kerneltype == "quadratic") {
    fx_quadratic <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      ED2 <- pdist::pdist(location[i, ],location[-i, ])@dist^2
      line_i[-i] <- 1 - ED2/(ED2 + (bandwidth))
      ind_i<- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)),
                  "val_i" = line_i[ind_i]))
    }

    results <- parallel::mclapply(1:dim(location)[1], fx_quadratic, mc.cores = ncores)
    tib <- tidyr::tibble(results)  %>%  tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]),
                                     x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                                    dim(location)[1]))
  }
  return(K_sparse)
}


#' @title Build a categorical/compositional kernel matrix.
#' @description This function constructs a kernel matrix from a feature matrix
#' (e.g., a cell-type indicator or cell-type proportion matrix) using either the
#' Hamming or Euclidean distance. The kernel between two samples \eqn{i} and \eqn{j}
#' is defined as \eqn{K_{ij} = \exp(-\gamma \cdot d(x_i, x_j))}, where
#' \eqn{d(\cdot,\cdot)} is either the Hamming distance (for discrete/categorical
#' features) or the squared Euclidean distance (for continuous features such as
#' cell-type proportions). Larger values of \eqn{\gamma} produce a kernel that
#' decays more rapidly with dissimilarity.
#' @param X An n by p feature matrix, where each row corresponds to a cell/spot
#' and each column to a feature.
#' @param gamma A positive numeric value controlling the decay rate of the kernel.
#' Default is 1.
#' @param disttype The type of distance to use in the kernel. Either \code{"hamming"}
#' (default; appropriate for categorical/indicator features) or \code{"euclidean"}
#' (appropriate for continuous features such as cell-type proportions).
#' @return An n by n symmetric positive semi-definite kernel matrix with diagonal
#' entries equal to 1.
Calculate_celltype_Kernel <- function(object, gamma = 1, disttype = "manhattan") {
  D <- as.matrix(stats::dist(as.matrix(object@cell_type_compositions), method = disttype))
  K <- exp(-gamma * D)
  object@cell_type_compositions_kernel<-K
  return(object)
}
