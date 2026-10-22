#' Rank-based inverse normal transform with zero imputation
#'
#' Zero counts are replaced by independent Uniform(0, 0.999) draws so that the
#' many tied zeros are broken at random. The values are then converted to normal
#' scores, \eqn{\Phi^{-1}((rank - 0.5) / n)}, and the scores are averaged over
#' `n_imp` imputations.
#'
#' @param y A numeric vector of counts for one gene.
#' @param n_imp Number of imputations to average over.
#' @return A numeric vector of the same length as `y`.
#' @examples
#' rank_int_impute(c(0, 0, 3, 1, 0, 7))
#' @export
rank_int_impute <- function(y, n_imp = 10) {
  n <- length(y)
  is_zero <- y == 0
  scores <- replicate(n_imp, {
    y_imp <- y
    y_imp[is_zero] <- stats::runif(sum(is_zero), 0, 0.999)
    stats::qnorm((rank(y_imp, ties.method = "average") - 0.5) / n)
  })
  rowMeans(matrix(scores, nrow = n))
}
