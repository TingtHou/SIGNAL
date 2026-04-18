#' @title Parallel KNN-based gene testing.
#' @description Applies rank-based normalization to each gene's expression vector
#' and tests each gene using MINQUE variance component estimation followed by a
#' chi-squared test on the specified components.
#' @param scaled_counts A gene-by-sample matrix of expression values.
#' @param KList_In A list of kernel matrices to be used for testing.
#' @param ncore Number of CPU cores to use for parallel processing.
#' @param TestingID Integer vector indicating which variance components to test.
#' @return A matrix with one row per gene; first column is the gene index,
#' remaining columns are the p-values for each component in \code{TestingID}.
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats qnorm
#' @importFrom KernelNN Linear MINQUE0 MNQTest0_Chi
#' @export
KNN_Parallel <- function(scaled_counts, KList_In, ncore, TestingID = c()) {

  # Build the kernel list used for variance component estimation
  HiddenKernel <- KernelNN::Linear(KList_In)

  results <- pbmcapply::pbmclapply(
    1:nrow(scaled_counts),
    function(i) {
      tryCatch({
        y <- scaled_counts[i, ]
        y <- stats::qnorm((rank(y, ties.method = "average") - 0.5) / length(y))

        res <- Testing_Each_Gene(y, HiddenKernel, i, TestingID)
        return(c(i, res))
      },
      error = function(e) {
        message("ERROR in gene ", i, ": ", e$message)
        return(c(i, rep(NA, length(TestingID) + 1)))
      })
    },
    mc.cores = ncore,
    mc.preschedule = FALSE,
    ignore.interactive = getOption("ignore.interactive", FALSE)
  )

  pMatrix <- do.call(rbind, results)
  return(pMatrix)
}

KNN<- function(scaled_counts, KList_In, ncore, TestingID = c()) {

  # Build the kernel list used for variance component estimation
  HiddenKernel <- KernelNN::Linear(KList_In)
  for (i in 1:nrow(scaled_counts)) {
    y <- scaled_counts[i, ]
    y <- stats::qnorm((rank(y, ties.method = "average") - 0.5) / length(y))

    res <- Testing_Each_Gene(y, HiddenKernel, i, TestingID)
  }


  return(pMatrix)
}


#' @title Test a single gene using MINQUE variance components.
#' @param y A normalized expression vector for one gene.
#' @param HiddenKernel A kernel list used for MINQUE estimation.
#' @param i Gene index (for error reporting).
#' @param TestingID Integer vector indicating which variance components to test.
#' @return A numeric vector of p-values on success, or NA on failure.
#'
#' @importFrom KernelNN MINQUE0 MNQTest0_Chi
Testing_Each_Gene <- function(y, HiddenKernel, i, TestingID) {
  tryCatch({
    vcs.result <- KernelNN::MINQUE0(KList = HiddenKernel, y = y)
    VCs        <- vcs.result$vcs
    Pvalue     <- KernelNN::MNQTest0_Chi(y, KList = HiddenKernel,
                                         vcs = VCs,
                                         ComponentID = TestingID)

    return(c(Pvalue$overall,Pvalue$components))
  },
  error = function(e) {
    cat(sprintf("Error at row %d: %s\n", i, conditionMessage(e)))
    return(rep(NA, length(TestingID)))
  })
}
