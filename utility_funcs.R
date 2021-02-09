require(Matrix)
require(matrixcalc)
library(mvtnorm)
require(rTensor)
require(abind)
require(foreach)

AR_model = function(m, rho) {
  #' @description Generate sparse inverse covariance matrix that has
  #' auto-regressive structure and its corresponding covariance matrix
  #' 
  #' @param m integer. the desired dimension of the cov matrix
  #' @param rho float. the weight parameter between 0 and 1 
  #'
  #' @return list. m x m inverse covariance matrix (inv) and its 
  #' corresponding covariance matrix (cov)
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname AR_model
  #' @export
  
  cov = matrix(0, nrow = m, ncol = m)
  inv = matrix(0, nrow = m, ncol = m)
  
  cov = rho^abs(row(cov) - col(cov))
  
  diag(inv) = 1+rho^2
  inv[1,1] = 1
  inv[m,m] = 1
  diag(inv[-1,]) = -rho
  diag(inv[,-1]) = -rho
  inv = inv / (1-rho^2)
  
  return(list(cov=cov, inv=inv))
}

SB_model = function(m, rho, num_subgraph=4, tol=1e-5){
  #' @importFrom Matrix
  #' @description Generate sparse inverse covariance matrix that has
  #' star-block structure and its corresponding covariance matrix
  #' 
  #' @param m integer. the desired dimension of the cov matrix
  #' @param rho float. the weight parameter between 0 and 1 
  #' @param num_subgraph integer. the desired number of star structures
  #' @param tol float. thresholding value for zeros
  #'
  #' @return list. m x m inverse covariance matrix (inv) and its 
  #' corresponding covariance matrix (cov)
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname SB_model
  #' @export
  
  size_subgraph = floor(m/num_subgraph)
  unequal_graph = (m != size_subgraph * num_subgraph)
  
  list_subgraph = list()
  list_invsubgraph = list()
  for (i in 1:num_subgraph){
    if(i==num_subgraph && unequal_graph){
      size_subgraph = m - size_subgraph*(i-1)
    }
    central_node = sample(1:size_subgraph,1)
    subgraph = matrix(rho^2, nrow=size_subgraph, ncol = size_subgraph)
    subgraph[,central_node] = rho
    subgraph[central_node,] = rho
    diag(subgraph)=1
    inv_subgraph = list(solve(subgraph))
    subgraph = list(subgraph)
    list_subgraph = c(list_subgraph, subgraph)
    list_invsubgraph = c(list_invsubgraph, inv_subgraph)
  }
  
  B = as.matrix(bdiag(list_subgraph))
  B[abs(B)<tol] = 0
  B_inv = as.matrix(bdiag(list_invsubgraph))
  B_inv[abs(B_inv)<tol] = 0
  
  return(list(cov=B, inv=B_inv))
}

ER_model = function(m, dens, wmin = 0.1, wmax = 0.3, corr=TRUE, tol=1e-5, tau_B=NA){
  #' @description Generate sparse inverse covariance matrix that has
  #' Erdos-Renyi model and its corresponding covariance matrix
  #' 
  #' @param m integer. the desired dimension of the cov matrix
  #' @param dens integer. the desired number of nonzero edges
  #' @param wmin float. the minimum weight between 0 and 1 for the edges 
  #' @param wmax float. the maximum weight between 0 and 1 for the edges 
  #' @param corr bool. TRUE for returning a correlation matrix
  #' @param tol float. thresholding value for zeros
  #'
  #' @return list. m x m inverse covariance matrix (inv) and its 
  #' corresponding covariance matrix (cov)
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname ER_model
  #' @export
  
  B_inv = 0.25*diag(m)
  row_sample = sample(1:m, size = dens)
  col_sample = sample(1:m, size = dens)
  rand_weight = runif(dens, min = wmin, max = wmax)
  
  for (i in 1:dens) {
    B_inv[row_sample[i], col_sample[i]] = B_inv[row_sample[i], col_sample[i]] - rand_weight[i]
    B_inv[col_sample[i], row_sample[i]] = B_inv[row_sample[i], col_sample[i]]
    B_inv[row_sample[i], row_sample[i]] = B_inv[row_sample[i], row_sample[i]] + rand_weight[i]
    B_inv[col_sample[i], col_sample[i]] = B_inv[col_sample[i], col_sample[i]] + rand_weight[i]
  }
  
  B = solve(B_inv)
  
  
  if(corr){
    diag_B = 1/sqrt(diag(B))
    B = diag(diag_B) %*% B %*% diag(diag_B)
    B_inv = solve(B)
  }
  
  B[abs(B)<tol] = 0
  B_inv[abs(B_inv)<tol] = 0
  return(list(cov=B, inv=B_inv))
}

##
Unif_cov_model = function(dim, min_eigen=1e-3){
#' @description Creates a covariance matrix with uniformly distributed weights
#' 
#' @param dim integer. The desired dimension of the covariance matrix
#' @param min_eigen float. The minimum bound for the eigenvalues of the 
#' desired covariance matrix
#'
#' @return dim x dim covariance matrix with uniformly distributed weights
#' @author Byoungwook Jang and Wayne Wang
#' @rdname generate_cov
#' @export

    mat_unif = matrix(runif(dim^2,-1,1),dim,dim)
    cov = t(mat_unif) %*% mat_unif
    cov = cov + diag(max(-1.2*min(eigen(cov)$values), min_eigen), dim)
    return(cov)
}

# generate sparse covariance
Unif_model = function(dim, dens, min_eigen = 1e-3){
#' @description Generate sparse inverse covariance matrix with 
#' uniformly distributed weights
#' 
#' @param dim integer. The desired dimension of the covariance matrix
#' @param dens integer. The sparsity level of the final inverse covariance
#' @param min_eigen float. The minimum bound for the eigenvalues of the 
#' desired covariance matrix
#'
#' @return dim x dim covariance matrix with dens number of entries
#' @author Byoungwook Jang and Wayne Wang
#' @rdname generate_inv_cov
#' @export

    U = as.matrix(rsparsematrix(dim,dim,dens))
    for (i in 1:nrow(U)) {
        for (j in 1:ncol(U)) {
        if (U[i,j]!=0){
            U[i,j] = sample(c(-1,1),1)
        }
        }
    }
    U = t(U)%*%U
    d_U = diag(U)
    U = apply(apply(U-diag(d_U), 1:2, function(x) min(x,1)), 1:2, function(x) max(x,-1))
    A = U + diag(d_U+1)
    A = A + diag(max(-1.2*min(eigen(A)$values),0.001),dim,dim)
    return(A)
}

# Calculates FNR, FPR, and
metric_edge = function(true_edge, est_edge){
  
  true_edge = true_edge*1
  est_edge = est_edge*1
  
  lower_tri = lower.tri(true_edge, diag = FALSE)
  pos = true_edge[lower_tri] != 0
  
  FP = sum((est_edge[lower_tri][!pos] != 0))
  FN = sum((est_edge[lower_tri][pos] == 0))
  TP = sum((est_edge[lower_tri][pos] != 0))
  TN = sum((est_edge[lower_tri][!pos] == 0))
  
  FPR = FP / (FP + TN)
  FNR = FN / (TP + FN)
  Precision = TP / (TP + FP)
  Recall = TP / (TP + FN)
  MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  metric_result = c(FP, FN, TP, TN, FPR, FNR, Precision, Recall, MCC)
  names(metric_result) = c('FP', 'FN', 'TP', 'TN', 'FPR', 'FNR', 'Precision', 'Recall', 'MCC')
  
  return(metric_result)
}

# Kronecker operations
kway_kronecker_sum <- function(matrices_list){
  #' @importFrom 
  #' @description K-way Kronecker sum
  #' @param matrices_list a list of matrices 
  #' @return dims x dims matrix
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname kway_kronecker_sum
  #' @export
  K <- length(matrices_list)
  if (K == 1) {
    return(matrices_list[[1]])
  } else{
    A <- kronecker_sum(matrices_list[[1]],matrices_list[[2]])
    if (K == 2) {
      return(A)
    } else {
      for (k in 3:K){
        A <- kronecker_sum(A,matrices_list[[k]])
      }
      return(A)
    }
  }
}

kronecker_sum_sparse = function(A,B){
  m = NCOL(A)
  n = NCOL(B)
  return(as(kronecker(Diagonal(n),A) + kronecker(B,Diagonal(m)), "dgCMatrix"))
}

kway_kronecker_sum_sparse = function(matrices_list){
  K <- length(matrices_list)
  if (K == 1) {
    return(matrices_list[[1]])
  } else{
    A <- kronecker_sum_sparse(matrices_list[[1]],matrices_list[[2]])
    if (K == 2) {
      return(A)
    } else {
      for (k in 3:K){
        A <- kronecker_sum_sparse(A,matrices_list[[k]])
      }
      return(A)
    }
  }
}

kway_kronecker = function(matrices_list){
  K <- length(matrices_list)
  if (K == 1) {
    return(matrices_list[[1]])
  } else{
    A <- kronecker(matrices_list[[1]],matrices_list[[2]])
    if (K == 2) {
      return(A)
    } else {
      for (k in 3:K){
        A <- kronecker(A,matrices_list[[k]])
      }
      return(A)
    }
  }
}

kmode_prod <- function(X,A,k){
  #' @importFrom 
  #' @description # k-mode product of a tensor and a matrix
  #' 
  #' @param X tensor of dims m_1,...,m_K
  #' @param A matrix of dimention J x m_k
  #' @param k mode
  #' @return a tensor of dims m_1,..m_k-1,J,m_k+1,...m_K
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname kmode_prod
  #' @export
  X_k_A <- k_fold(mat_mult(A,k_unfold(X,m=k)@data),m=k,modes=X@modes)
  return(X_k_A)
}

# function to update the weights
update_diag = function(N,dims,X_N,Y_N){
  #' @importFrom 
  #' @description # k-mode product of a tensor and a matrix
  #' 
  #' @param N sample size
  #' @param dims m_1,...m_K
  #' @param X_N last mode unfolding of the data tensor
  #' @param Y_N a tensor of the same dims of the data
  #' @return a tensor contains all the weights
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname update_diag
  #' @export
  a = (1/N)*diag(mat_mult(t(X_N),X_N))
  b = (1/N)*diag(mat_mult(t(X_N),Y_N))
  w = (-b+sqrt(b^2+4*a))/(2*a)
  return(w)
}

# generate sylvester data
gen_sylvester_data = function(N, dims, Omega_sqrt){
  #' @importFrom
  #' @description Generate sylvester matrix (standardized) normal variable with relationship
  #' AX + XB =  C
  #'
  #' @param N number of replicates
  #' @param dims tensor dimensions
  #' @param Omega_sqrt prod(dims) x prod(dims) square root inverse covariance matrix
  #' @return dims x N random tensor corresponds to N samples of dims matrix
  #' @author Byoungwook Jang and Wayne Wang
  #' @rdname gen_sylvester_data
  #' @export

  Omega_sqrt <- as(Omega_sqrt, "dgCMatrix")

  if (N == 1){
    Z = t(stdmvrnormArma(N,prod(dims)))
    X = mvrnormEigen(Z,Omega_sqrt)
    X = array(X[,1],dim = c(dims,1))
    X = as.tensor(X)
  } else{
    Z = t(stdmvrnormArma(N,prod(dims)))
    X = mvrnormEigen(Z,Omega_sqrt)
    X = apply(X,1,scale)
    X = foreach(i=1:nrow(X),.combine=function(x,y) abind(x,y,rev.along=1)) %do% array(X[i,],dim=c(dims,1))
    X = as.tensor(X)
  }

  return(X)
}

# convert a inv covariance matrix to a partial correlation matrix
invcov2parcor <- function(invcov){
  p <- ncol(invcov)
  parcor <- diag(p)
  
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      parcor[i,j] <- -invcov[i,j]/sqrt(invcov[i,i]*invcov[j,j])
      parcor[j,i] <- parcor[i,j]
    }
  }
  
  return(parcor)
}
