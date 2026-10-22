# Core variance component tests.
#
# For kernels K_1, ..., K_m (K_m = identity), MINQUE(0) estimates are
# tau = C^{-1} u with C_ij = tr(K_i K_j) and u_j = y' K_j y, so each estimate
# is a quadratic form tau_i = y' Q_i y with Q_i = sum_j (C^{-1})_ij K_j.
#
# * Overall test (tau_1 = ... = tau_{m-1} = 0): with r = sum(tau_effects) / tau_m,
#   p = P(z' (B - r D) z > 0), B = sum of Q_i over the effects, D = Q_m.
# * Component test for tau_k: p = P(y' Q_k y > tau_k) with y ~ N(0, V0), where
#   V0 = sum_{j != k} tau0_j K_j is the MINQUE(0) fit without kernel k.
#
# Both tail probabilities use Liu's approximation, which only needs the power
# sums tr(A^h), h = 1..4, of the matrix A of the quadratic form. These traces
# are polynomials in r (or tau0) whose coefficients are traces of kernel
# products, so they are computed once per kernel set and evaluated for all
# genes at once ("moment" method). The "eigen" method computes the eigenvalues
# of A for every gene and applies Liu's approximation to them. The "davies"
# method uses those eigenvalues with Davies' algorithm (CompQuadForm::davies)
# for every test. The default "davies_liu" method is the original SIGNAL code
# used for the real-data analyses: Davies for the overall test and Liu for the
# component tests.

# Liu, Tang and Zhang (2009) approximation to P(Q > q), vectorised, written in
# terms of the power sums c_h = sum(lambda^h). Matches CompQuadForm::liu().
.liu <- function(q, c1, c2, c3, c4) {
  s1 <- c3 / c2^1.5
  s2 <- c4 / c2^2
  tstar <- (q - c1) / sqrt(2 * c2)
  big <- !is.na(s1) & !is.na(s2) & s1^2 > s2
  a <- ifelse(big, 1 / (s1 - sqrt(pmax(s1^2 - s2, 0))), 1 / s1)
  delta <- ifelse(big, s1 * a^3 - a^2, 0)
  l <- ifelse(big, a^2 - 2 * delta, c2^3 / c3^2)
  suppressWarnings(
    stats::pchisq(tstar * sqrt(2) * a + l + delta, df = l, ncp = delta,
                  lower.tail = FALSE)
  )
}

.liu_eigen <- function(q, lambda) {
  .liu(q, sum(lambda), sum(lambda^2), sum(lambda^3), sum(lambda^4))
}

# Davies (1980) P(Q > q). Its accuracy is limited to about `acc`; values it
# returns below 0 are set to 0.
.davies_eigen <- function(q, lambda) {
  p <- CompQuadForm::davies(q, lambda)$Qq
  min(max(p, 0), 1)
}

# Eigenvalues of Q %*% V (Q, V symmetric). If V = R'R is positive definite they
# equal the eigenvalues of the symmetric matrix R Q R', which are much faster to
# compute; otherwise use the real parts of the eigenvalues of Q %*% V.
.product_eigenvalues <- function(Q, V) {
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (!is.null(R)) {
    return(eigen(R %*% Q %*% t(R), symmetric = TRUE, only.values = TRUE)$values)
  }
  Re(eigen(Q %*% V, only.values = TRUE)$values)
}

# tr(X %*% Y) without forming the product.
.tr2 <- function(X, Y) sum(X * t(Y))

# u_j = y' K_j y for every gene (rows of Y), processed in chunks of genes.
.quad_forms <- function(Y, Klist, chunk_size = 1000) {
  G <- nrow(Y)
  U <- matrix(NA_real_, G, length(Klist))
  for (start in seq(1, G, by = chunk_size)) {
    idx <- start:min(G, start + chunk_size - 1)
    Yc <- Y[idx, , drop = FALSE]
    for (j in seq_along(Klist)) {
      U[idx, j] <- rowSums((Yc %*% Klist[[j]]) * Yc)
    }
  }
  U
}

.linear_comb <- function(Klist, w) {
  Reduce(`+`, Map(`*`, Klist, w))
}

# Power-sum coefficients for the overall test, A = B - r D (B, D symmetric).
.overall_moments <- function(B, D) {
  BB <- B %*% B
  BD <- B %*% D
  DD <- D %*% D
  list(
    B = sum(diag(B)), D = sum(diag(D)),
    BB = sum(B * B), BD = sum(B * D), DD = sum(D * D),
    BBB = sum(BB * B), BBD = sum(BB * D), BDD = sum(DD * B), DDD = sum(DD * D),
    BBBB = sum(BB * BB), BBBD = .tr2(BB, BD), BBDD = sum(BB * DD),
    BDBD = .tr2(BD, BD), BDDD = sum(BD * DD), DDDD = sum(DD * DD)
  )
}

.overall_pvalue_moment <- function(r, mo) {
  c1 <- mo$B - r * mo$D
  c2 <- mo$BB - 2 * r * mo$BD + r^2 * mo$DD
  c3 <- mo$BBB - 3 * r * mo$BBD + 3 * r^2 * mo$BDD - r^3 * mo$DDD
  c4 <- mo$BBBB - 4 * r * mo$BBBD + r^2 * (4 * mo$BBDD + 2 * mo$BDBD) -
    4 * r^3 * mo$BDDD + r^4 * mo$DDDD
  .liu(0, c1, c2, c3, c4)
}

# Trace tensors tr(M_j1 ... M_jh), h = 1..4, for M_j = Q_k K_j.
.component_moments <- function(Qk, Knull) {
  M <- lapply(Knull, function(K) Qk %*% K)
  J <- length(M)
  T1 <- vapply(M, function(X) sum(diag(X)), numeric(1))
  T2 <- matrix(0, J, J)
  T3 <- array(0, c(J, J, J))
  T4 <- array(0, c(J, J, J, J))
  P <- vector("list", J * J)
  for (j in seq_len(J)) for (l in seq_len(J)) {
    T2[j, l] <- .tr2(M[[j]], M[[l]])
    P[[(l - 1) * J + j]] <- M[[j]] %*% M[[l]]
  }
  Mt <- lapply(M, t)
  rm(M)
  for (j in seq_len(J)) for (l in seq_len(J)) {
    Pjl <- P[[(l - 1) * J + j]]
    for (m in seq_len(J)) {
      T3[j, l, m] <- sum(Pjl * Mt[[m]])
      for (p in seq_len(J)) {
        T4[j, l, m, p] <- .tr2(Pjl, P[[(p - 1) * J + m]])
      }
    }
  }
  list(T1 = T1, T2 = T2, T3 = T3, T4 = T4)
}

.component_pvalue_moment <- function(q, A, mo) {
  J <- ncol(A)
  c1 <- drop(A %*% mo$T1)
  c2 <- rowSums((A %*% mo$T2) * A)
  c3 <- c4 <- numeric(nrow(A))
  for (j in seq_len(J)) for (l in seq_len(J)) for (m in seq_len(J)) {
    a3 <- A[, j] * A[, l] * A[, m]
    c3 <- c3 + mo$T3[j, l, m] * a3
    for (p in seq_len(J)) {
      c4 <- c4 + mo$T4[j, l, m, p] * a3 * A[, p]
    }
  }
  .liu(q, c1, c2, c3, c4)
}

# Run the overall test and the component tests in `test_index` for all genes
# (rows of Y) with one kernel list (last kernel = identity).
# Returns a G x (1 + length(test_index)) matrix of p-values.
.signal_engine <- function(Y, Klist, test_index = 1:3,
                           method = c("davies_liu", "davies", "moment", "eigen")) {
  method <- match.arg(method)
  # Tail probability used with per-gene eigenvalues: overall test, component tests
  overall_tail <- if (method %in% c("davies_liu", "davies")) .davies_eigen else .liu_eigen
  component_tail <- if (method == "davies") .davies_eigen else .liu_eigen
  m <- length(Klist)
  G <- nrow(Y)
  out <- matrix(NA_real_, G, 1 + length(test_index))

  C <- matrix(0, m, m)
  for (i in seq_len(m)) for (j in i:m) {
    C[i, j] <- C[j, i] <- sum(Klist[[i]] * Klist[[j]])
  }
  Cinv <- tryCatch(solve(C), error = function(e) NULL)
  if (is.null(Cinv) || anyNA(Cinv) || any(diag(Cinv) < 0)) {
    warning("Kernel Gram matrix is singular or not positive definite; ",
            "returning NA for this kernel.", call. = FALSE)
    return(out)
  }
  Q <- lapply(seq_len(m), function(i) .linear_comb(Klist, Cinv[i, ]))
  U <- .quad_forms(Y, Klist)
  tau <- U %*% t(Cinv)

  effects <- seq_len(m - 1)
  r <- rowSums(tau[, effects, drop = FALSE]) / tau[, m]
  B <- Reduce(`+`, Q[effects])
  D <- Q[[m]]
  if (method == "moment") {
    out[, 1] <- .overall_pvalue_moment(r, .overall_moments(B, D))
  } else {
    out[, 1] <- vapply(seq_len(G), function(g) {
      lambda <- eigen(B - r[g] * D, symmetric = TRUE, only.values = TRUE)$values
      overall_tail(0, lambda)
    }, numeric(1))
  }
  rm(B, D)

  for (t in seq_along(test_index)) {
    k <- test_index[t]
    null_set <- setdiff(seq_len(m), k)
    Cnull_inv <- tryCatch(solve(C[null_set, null_set]), error = function(e) NULL)
    if (is.null(Cnull_inv)) next
    tau0 <- U[, null_set, drop = FALSE] %*% t(Cnull_inv)
    if (method %in% c("moment", "davies_liu")) {
      mo <- .component_moments(Q[[k]], Klist[null_set])
      p <- .component_pvalue_moment(tau[, k], tau0, mo)
      if (method == "davies_liu") {
        # Moments equal the power sums of the eigenvalues of Q_k V0 when V0 is
        # positive definite (all eigenvalues real). Otherwise use the real parts
        # of the eigenvalues, as the original code does.
        for (g in seq_len(G)) {
          V0 <- .linear_comb(Klist[null_set], tau0[g, ])
          if (is.null(tryCatch(chol(V0), error = function(e) NULL))) {
            lambda <- Re(eigen(Q[[k]] %*% V0, only.values = TRUE)$values)
            p[g] <- .liu_eigen(tau[g, k], lambda)
          }
        }
      }
      out[, t + 1] <- p
    } else {
      out[, t + 1] <- vapply(seq_len(G), function(g) {
        V0 <- .linear_comb(Klist[null_set], tau0[g, ])
        component_tail(tau[g, k], .product_eigenvalues(Q[[k]], V0))
      }, numeric(1))
    }
  }
  out
}

#' Cauchy combination of p-values
#'
#' Combines the p-values in each row (one column per spatial kernel) with the
#' Cauchy combination test (ACAT; Liu and Xie, 2020). This is the
#' `CombinePValues()` function of the original SIGNAL analysis code: p-values
#' of exactly 0 are set to 5.55e-17 and those within 1e-3 of 1 are set to 0.99
#' before combining, and combined p-values of 0 are set to 5.55e-17.
#'
#' @param pvalues A matrix of p-values, one row per gene.
#' @param weights Optional matrix of weights with the same dimensions as
#'   `pvalues`; by default every column gets weight 1 / ncol(pvalues).
#' @param na.rm If `FALSE` (default, as in the original code) a row with any
#'   missing p-value gives a missing combined p-value. If `TRUE`, missing
#'   p-values are dropped and the remaining weights are rescaled to sum to 1.
#' @return A numeric vector of combined p-values.
#' @examples
#' P <- matrix(c(0.01, 0.2, 0.5, 0.04, NA, 0.9), nrow = 2)
#' combine_pvalues(P)
#' combine_pvalues(P, na.rm = TRUE)
#' @export
combine_pvalues <- function(pvalues, weights = NULL, na.rm = FALSE) {
  if (!is.matrix(pvalues)) pvalues <- as.matrix(pvalues)
  pvalues[which(pvalues == 0)] <- 5.55e-17
  pvalues[which((1 - pvalues) < 1e-3)] <- 0.99

  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if (is.null(weights)) {
    weights <- matrix(rep(1.0 / num_pval, num_pval * num_gene), ncol = num_pval)
  }
  if (nrow(weights) != num_gene || ncol(weights) != num_pval) {
    stop("the dimensions of weights does not match that of combined pvalues")
  }

  Cstat <- tan((0.5 - pvalues) * pi)
  wCstat <- weights * Cstat
  if (na.rm) {
    w <- weights
    w[is.na(Cstat)] <- 0
    Cbar <- rowSums(wCstat, na.rm = TRUE) / rowSums(w)
  } else {
    Cbar <- apply(wCstat, 1, sum)
  }
  combined_pval <- 1.0 - stats::pcauchy(Cbar)
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
  combined_pval
}
