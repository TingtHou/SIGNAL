# Reference implementation copied from the original analysis code
# (SingleCell/KNN/MINQUE.R, KNN/MNQTestChi.R, KNN/KNN.R), used to check that
# the package reproduces it. Requires CompQuadForm.

ref_MINQUE0 <- function(KList, y) {
  VCs <- length(KList)
  C.mat <- matrix(NA, nrow = VCs, ncol = VCs)
  RightY <- c()
  for (i in 1:VCs) {
    for (j in i:VCs) C.mat[j, i] <- C.mat[i, j] <- sum(t(KList[[i]]) * KList[[j]])
    RightY <- c(RightY, t(y) %*% KList[[i]] %*% y)
  }
  list(vcs = solve(C.mat) %*% RightY)
}

ref_Q <- function(KList) {
  VCs <- length(KList)
  C.mat <- matrix(NA, nrow = VCs, ncol = VCs)
  for (i in 1:VCs) for (j in i:VCs) C.mat[j, i] <- C.mat[i, j] <- sum(t(KList[[i]]) * KList[[j]])
  C.mat <- solve(C.mat)
  lapply(1:VCs, function(i) Reduce("+", mapply("*", KList, C.mat[i, ], SIMPLIFY = FALSE)))
}

ref_overall <- function(KList, est_vcs, index_interest, pfun = CompQuadForm::liu) {
  VCs <- length(est_vcs)
  QList <- ref_Q(KList)
  ratio <- sum(est_vcs[index_interest]) / est_vcs[VCs]
  sigma_total <- Reduce("+", QList[index_interest]) - ratio * QList[[VCs]]
  eigen_test <- eigen(sigma_total, symmetric = TRUE, only.values = TRUE)$values
  pfun(q = 0, lambda = eigen_test)
}

ref_single <- function(KList, est_vcs, est_vcs_h0, ParaTestingIndex, pfun = CompQuadForm::liu) {
  QList <- ref_Q(KList)
  ratio <- sum(est_vcs[ParaTestingIndex])
  sigma_total <- Reduce("+", QList[ParaTestingIndex])
  V0 <- Reduce("+", mapply("*", KList[-ParaTestingIndex], est_vcs_h0, SIMPLIFY = FALSE))
  eigen_test <- Re(eigen(sigma_total %*% V0, only.values = TRUE)$values)
  pfun(q = ratio, lambda = eigen_test)
}

# Testing_Each_Gene() with TestingID = list(1, 2, 3), as in SimluationImp.R.
ref_test_gene <- function(y, KList, pfun = CompQuadForm::liu, pfun_overall = pfun) {
  VCs <- ref_MINQUE0(KList, y)$vcs
  p <- ref_overall(KList, VCs, 1:(length(KList) - 1), pfun_overall)
  for (k in 1:3) {
    h0 <- ref_MINQUE0(KList[-k], y)$vcs
    p <- c(p, ref_single(KList, VCs, h0, k, pfun))
  }
  p
}

# Kernel construction exactly as in SimluationImp.R lines 270-328.
ref_kernels <- function(KCellType, Kspatial) {
  KList_In <- list(KCellType, Kspatial,
                   (KCellType - mean(KCellType)) * (Kspatial - mean(Kspatial)))
  KList_In <- lapply(KList_In, function(K) {
    matrix(scale(as.vector(K), scale = FALSE), nrow = nrow(K), ncol = ncol(K))
  })
  c(KList_In, list(diag(nrow(Kspatial))))
}

# CombinePValues() verbatim from SingleCell/funcs/utilities.R.
CombinePValues <- function(pvalues, weights=NULL){
  if(!is.matrix(pvalues)){pvalues <- as.matrix(pvalues)}
  ## to avoid extremely values
  pvalues[which(pvalues==0)] <- 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] <- 0.99
  
  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if(is.null(weights)){
    weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }# end fi
  if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
    stop("the dimensions of weights does not match that of combined pvalues")
  }# end fi
  
  Cstat <- tan((0.5 - pvalues)*pi)
  
  wCstat <- weights*Cstat
  Cbar <- apply(wCstat, 1, sum)
  #combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval <- 1.0 - pcauchy(Cbar)	
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
  return(combined_pval)
}# end func
