library(rTensor)
library(matrixcalc)
library(Matrix)
library(abind)
library(foreach)
library(Rcpp)
sourceCpp("helper_funcs.cpp")
source("utility_funcs.R")


syglasso <- function(X,lambda,Psi_true=NA,Omega_true=NA,maxiter=100,epsilon=1e-6){
  # assuming X is in tensor format
  K <- X@num_modes - 1
  dims <- X@modes[1:K]
  
  # pre-compute all mode-k gram matrices and fibers
  X_N <- k_unfold(X,m=K+1)@data 
  Grams <- list()
  X_fibers <- list()
  for (k in 1:K){
    X_fibers <- append(X_fibers,list(k_unfold(X,m=k)@data))
    Grams <- append(Grams,list(mat_mult(X_fibers[[k]],t(X_fibers[[k]]))))
  }
  
  # initialization 
  Psi_init <- list()
  for (k in 1:K){
    # identity init
    Psi_init <- append(Psi_init,list(diag(dims[k])))
  }
  w_init <- as.tensor(array(diag(kway_kronecker_sum(lapply(1:K, function(k) diag(diag(Psi_init[[k]]))))),dim=dims)) # took long
  Omega_init <- mat_mult(kway_kronecker_sum(Psi_init)+diag(rTensor::vec(w_init)),t(kway_kronecker_sum(Psi_init)+diag(rTensor::vec(w_init))))
 
  dist_iter_opt_A = rep(NA, maxiter)
  dist_iter_true_A = rep(NA, maxiter)
  dist_iter_opt_B = rep(NA, maxiter)
  dist_iter_true_B = rep(NA, maxiter)  
  dist_iter_opt_Omega = rep(NA, maxiter)
  dist_iter_true_Omega = rep(NA, maxiter)
  
  Psi_current_seq <- list()
  Omega_current_seq <- list()
  
  time_seq <- numeric()
  
  Psi_current <- Psi_init
  w_current <- rTensor::vec(w_init)
  Omega_current <- Omega_init
  iter <- 1
  converged <- F
  
  while ((!converged)&(iter<maxiter)){
    Omega_old <- Omega_current
    Psi_old <- Psi_current
    W <- array(w_current,dim=c(dims,1))
    W <- foreach(i=1:N,.combine=function(x,y) abind(x,y,rev.along=1)) %do% W
    W <- as.tensor(W)
    
    # update off-diagonal elements for each matrices
    for (k in 1:K){
      Psi_k <- Psi_current[[k]]
      m_k <- ncol(Psi_k)
      X_k <- X_fibers[[k]]
      gram_k <- Grams[[k]]
      W_k <- k_unfold(W,m=k)@data
      W_gram <- mat_mult(hadamard.prod(W_k,X_k),t(X_k))
      
      # pre-compute quantities needed for updating eqn
      sum_gram_othermodes <- matrix(0,m_k,m_k)
      for (l in 1:K){
        if (l == k){next}
        sum_gram_othermodes <- sum_gram_othermodes + 
          mat_mult(X_k,t(k_unfold(kmode_prod(X,Psi_current[[l]]-diag(diag(Psi_current[[l]])),l),m=k)@data))
      }
      
      # updating eqn
      Psi_k <- update_off_diag(N,m_k,lambda,Psi_k,gram_k,X_k,sum_gram_othermodes,W_gram)
      Psi_current[[k]] <- Psi_k
    }
    
    # update weights(diagonal elements) 
    # pre-compute quantities needed
    sum_kmodes_prod <- array(0,dim=c(dims,N))
    for (l in 1:K){
      sum_kmodes_prod <- sum_kmodes_prod +
        kmode_prod(X,Psi_current[[l]]-diag(diag(Psi_current[[l]])),l)@data
    }
    sum_kmodes_prod <- as.tensor(sum_kmodes_prod)
    sum_kmodes_prod_N <- k_unfold(sum_kmodes_prod,m=K+1)@data
    
    # updating eqn
    w_current <- update_diag(N,dims,X_N,sum_kmodes_prod_N)
    
    # update Omega (if not too large in size)
    # Omega_current <- mat_mult(kway_kronecker_sum(Psi_current)+diag(w_current),t(kway_kronecker_sum(Psi_current)+diag(w_current)))
    # dist_iter_Omega[iter] <- norm(Omega_current-Omega_old, type='m')
    # dist_true_Omega[iter] <- norm(Omega_true-Omega_current, type='f')
    # Omega_current_seq <- append(Omega_current_seq, list(Omega_current))
    # 
    # dist_true_A[iter] <- norm(Psi_true[[1]]-diag(diag(Psi_true[[1]]))-Psi_current[[1]], type='f') 
    # dist_true_B[iter] <- norm(Psi_true[[2]]-diag(diag(Psi_true[[2]]))-Psi_current[[2]], type='f')
    # 
    # dist_iter_A[iter] <- norm(Psi_current[[1]]-Psis_hat[[1]], type='f')/norm(Psis[[1]], type='f')
    # dist_iter_B[iter] <- norm(Psi_current[[2]]-Psis_hat[[2]], type='f')/norm(Psis[[2]], type='f')
    
    # update Psis
    Psi_current_seq <- append(Psi_current_seq, list(Psi_current))
    
    time_seq <- c(time_seq,proc.time()[3])
    
    # check convergence
    if (max(unlist(lapply(1:K, function(k) norm(Psi_current[[k]]-Psi_old[[k]],type='f'))))<epsilon){
      converged <- T
    } else{
      iter <- iter + 1
    }
    
    # remove unnecessary objects
    # rm(list=c('sum_kmodes_prod','sum_kmodes_prod_N','X_k','gram_k','W_k','W','sum_gram_othermodes'))
    gc()
  }
  # Omega_current <- mat_mult(kway_kronecker_sum(Psi_current)+diag(w_current),t(kway_kronecker_sum(Psi_current)+diag(w_current)))
  # optimization error
  # dist_iter_opt_A <- unlist(lapply(1:(iter-1), function(t) norm(Psi_current_seq[[t]][[1]] - Psi_current[[1]], type='f') / norm(Psi_true[[1]], type='f')))
  # dist_iter_opt_B <- unlist(lapply(1:(iter-1), function(t) norm(Psi_current_seq[[t]][[2]] - Psi_current[[2]], type='f') / norm(Psi_true[[2]], type='f')))
  # dist_iter_opt_Omega <- unlist(lapply(1:(iter-1), function(t) norm(Omega_current_seq[[t]] - Omega_current, type = 'f') / norm(Omega_true, type='f')))
  # statistical error
  # dist_iter_true_A <- unlist(lapply(1:(iter-1), function(t) norm(Psi_current_seq[[t]][[1]] - Psi_true[[1]] + diag(diag(Psi_true[[1]])), type='f') / norm(Psi_true[[1]], type='f')))
  # dist_iter_true_B <- unlist(lapply(1:(iter-1), function(t) norm(Psi_current_seq[[t]][[2]] - Psi_true[[2]] + diag(diag(Psi_true[[2]])), type='f') / norm(Psi_true[[2]], type='f')))
  # dist_iter_true_Omega <- unlist(lapply(1:(iter-1), function(t) norm(Omega_current_seq[[t]] - Omega_true, type = 'f') / norm(Omega_true, type='f')))
  
  return(list(Omega=Omega_current,Psi=Psi_current,W=w_current,
              Convergence=converged,Iter=iter,time_seq=time_seq,
              Psi_seq=Psi_current_seq,Omega_seq=Omega_current_seq,
              dist_iter_true_A=dist_iter_true_A,dist_iter_true_B=dist_iter_true_B,
              dist_iter_opt_A=dist_iter_opt_A,dist_iter_opt_B=dist_iter_opt_B,
              dist_iter_true_Omega=dist_iter_true_Omega,dist_iter_opt_Omega=dist_iter_opt_Omega))
}


