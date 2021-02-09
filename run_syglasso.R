rm(list=ls())
source("syglasso.R")

library(matrixcalc)
library(Matrix)
library(rTensor)
library(abind)
library(foreach)

# specify dimensions
dims<- c(10,10,10)
N <- 100
K <- length(dims)
dens <- 0.1

# create true sylvester parameters of types: AR, SB, ER, or Unif defined
# in the utility_funcs.R
Psis <- list()
for (k in 1:length(dims)){
  Psis <- append(Psis,list(AR_model(dims[k],0.25)$inv))
}

# simulate data
Omega_sqrt <- kway_kronecker_sum_sparse(Psis)
X <- gen_sylvester_data(N,dims,Omega_sqrt)

# run SyGlasso 
lambda <- 2.5*K*N*qnorm(1-0.05/(2*(prod(dims))^2))/sqrt(N) # penalty parameter
tmp <- proc.time()
result_tnsr <- syglasso(X,lambda,Omega_true=NA,Psi_true=NA,maxiter=50,epsilon=1e-6)
proc.time() - tmp

Psis_hat <- result_tnsr$Psi

Psis_hat_seq <- result_tnsr$Psi_seq

####################################################################################
# Check graph recovery performance
# Omega_hat <- mat_mult(kway_kronecker_sum(Psis_hat),t(kway_kronecker_sum(Psis_hat)))
# metric_edge(Omega,Omega_hat)
metric_edge(Psis[[1]],Psis_hat[[1]])
metric_edge(Psis[[2]],Psis_hat[[2]])
metric_edge(Psis[[3]],Psis_hat[[3]])
sum_metric <- metric_edge(Psis[[1]],Psis_hat[[1]])+metric_edge(Psis[[2]],Psis_hat[[2]])+metric_edge(Psis[[3]],Psis_hat[[3]])

#######################################################################################
# Check run time vs. performance
# Time vs. MCC
times_syglasso <- as.numeric(result_tnsr$time_seq - result_tnsr$time_seq[1])
avg_mcc_seq <- function(Psis_hat_seq,Psis){
  MCC <- numeric()
  for (i in 1:length(Psis_hat_seq)){
    Psis_hat <- Psis_hat_seq[[i]]
    sum_metric <- metric_edge(Psis[[1]],Psis_hat[[1]])+
      metric_edge(Psis[[2]],Psis_hat[[2]])+
      metric_edge(Psis[[3]],Psis_hat[[3]])
    TP <- sum_metric['TP']
    FP <- sum_metric['FP']
    TN <- sum_metric['TN']
    FN <- sum_metric['FN']
    MCC <- c(MCC,as.numeric((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  }
  return(MCC)
}
mcc_syglasso <- avg_mcc_seq(Psis_hat_seq,Psis)
plot(times_syglasso,mcc_syglasso,type = 'o')

# Time vs. NRMSE
nrmse_seq <- function(Psis_hat_seq,Psis){
  nrmse <- numeric()
  for (i in 1:length(Psis_hat_seq)){
    Psis_hat <- Psis_hat_seq[[i]]
    tmp <- 0
    for (k in 1:K){
      tmp <- tmp + norm(Psis[[k]]-diag(diag(Psis[[k]]))-Psis_hat[[k]], type = "f")/norm(Psis[[k]]-diag(diag(Psis[[k]])), type = "f")
    }
    nrmse <- c(nrmse,tmp/K)
  }
  return(nrmse)
}
nrmse_syglasso <- nrmse_seq(Psis_hat_seq,Psis)
plot(times_syglasso[-1],nrmse_syglasso[-1],type = 'o')

# Time vs MCC&NRMSE
plot(times_syglasso,avg_mcc,type = 'o',col='blue',ylim = c(0.3,1.5))
lines(times_syglasso[-1],nrmse_syglasso[-1],type = 'o',col='red')


########################################################################################
# Check statistical vs. optimization error
# F-norm for statistical error
plot(log(result_tnsr$dist_iter_true_A),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")
plot(log(result_tnsr$dist_iter_true_B),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")
plot(log(result_tnsr$dist_iter_true_Omega),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")

# F-norm for optimization error
plot(log(result_tnsr$dist_iter_opt_A),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")
plot(log(result_tnsr$dist_iter_opt_B),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")
plot(log(result_tnsr$dist_iter_opt_Omega),type='o',xlab="Iterations",ylab="Relative Frobenius Norm")
