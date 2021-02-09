// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath> 
#include <algorithm>
#include <Eigen/SparseCholesky>

using namespace Rcpp;

// // [[Rcpp::export]]
// Eigen::MatrixXd update_off_diag(int N, int m_k, double lambda, Eigen::MatrixXd Psi_k, Eigen::MatrixXd gram_k, Eigen::MatrixXd X_k, Eigen::MatrixXd gram_l) {
// 
//   for (int p = 0; p < m_k - 1; p++){
//     for (int q = p + 1; q < m_k; q++){
//       Psi_k(p, q) = 0;
//       Psi_k(q, p) = 0;
//       Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k;
//       double num = (-1.0/N) * (gram_Psi_k(q, p) + gram_Psi_k(p, q) + gram_l(q, p) + gram_l(p, q));
//       double denom = (1.0/N) * (gram_k(q, q) + gram_k(p, p));
//       Psi_k(p, q) = (num/std::abs(num)) * std::max(std::abs(num) - lambda/N, 0.0)/denom;
//       Psi_k(q, p) = Psi_k(p, q);
//     }
//   }
// 
//   return Psi_k;
// }
// 
// // [[Rcpp::export]]
// Eigen::MatrixXd update_diag(int N, int m_k, int d, double sum_trace, Eigen::MatrixXd Psi_k, Eigen::MatrixXd gram_k, Eigen::MatrixXd X_k, Eigen::MatrixXd gram_l) {
// 
//   for (int p = 0; p < m_k; p++){
//     Psi_k(p, p) = 0;
//     Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k;
//     double a = (1.0/N) * gram_k(p, p);
//     double b = (1.0/N) * (gram_k(p, p) * sum_trace/(d/m_k) + gram_Psi_k(p, p) + gram_l(p, p));
//     double c = (1.0/N) * ((gram_Psi_k(p, p) + gram_l(p, p))*sum_trace/(d/m_k) - d/m_k);
//     Psi_k(p, p) = (-b + sqrt(pow(b, 2.0) - 4.0*a*c))/(2.0*a);
//   }
// 
//   return Psi_k;
// }

// // [[Rcpp::export]]
// Eigen::MatrixXd update_off_diag(int N, int m_k, double lambda, Eigen::Map<Eigen::MatrixXd> Psi_k, const Eigen::Map<Eigen::MatrixXd> gram_k, Eigen::Map<Eigen::MatrixXd> X_k, Eigen::Map<Eigen::MatrixXd> gram_l) {
// 
//   for (int p = 0; p < m_k - 1; p++){
//     for (int q = p + 1; q < m_k; q++){
//       Psi_k(p, q) = 0;
//       Psi_k(q, p) = 0;
//       Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k;
//       double num = (-1.0/N) * (gram_Psi_k(q, p) + gram_Psi_k(p, q) + gram_l(q, p) + gram_l(p, q));
//       double denom = (1.0/N) * (gram_k(q, q) + gram_k(p, p));
//       Psi_k(p, q) = (num/std::abs(num)) * std::max(std::abs(num) - lambda/N, 0.0)/denom;
//       Psi_k(q, p) = Psi_k(p, q);
//     }
//   }
// 
//   return Psi_k;
// }
// 
// // [[Rcpp::export]]
// Eigen::MatrixXd update_diag(int N, int m_k, int d, double sum_trace, Eigen::Map<Eigen::MatrixXd> Psi_k, const Eigen::Map<Eigen::MatrixXd> gram_k, Eigen::Map<Eigen::MatrixXd> X_k, Eigen::Map<Eigen::MatrixXd> gram_l) {
// 
//   for (int p = 0; p < m_k; p++){
//     Psi_k(p, p) = 0;
//     Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k;
//     double a = (1.0/N) * gram_k(p, p);
//     double b = (1.0/N) * (gram_k(p, p) * sum_trace/(d/m_k) + gram_Psi_k(p, p) + gram_l(p, p));
//     double c = (1.0/N) * ((gram_Psi_k(p, p) + gram_l(p, p))*sum_trace/(d/m_k) - d/m_k);
//     Psi_k(p, p) = (-b + sqrt(pow(b, 2.0) - 4.0*a*c))/(2.0*a);
//   }
// 
//   return Psi_k;
// }

// // [[Rcpp::export]]
// Eigen::MatrixXd update_diag(int N, int m_k, int d, double sum_trace, Eigen::Map<Eigen::MatrixXd> Psi_k, const Eigen::Map<Eigen::MatrixXd> gram_k, Eigen::Map<Eigen::MatrixXd> X_k, Eigen::Map<Eigen::MatrixXd> gram_l) {
//   
//   for (int p = 0; p < m_k; p++){
//     Psi_k(p, p) = 0;
//     Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k;
//     double a = (1.0/N) * gram_k(p, p);
//     double b = (1.0/N) * (gram_k(p, p) * sum_trace/(d/m_k) + gram_Psi_k(p, p) + gram_l(p, p));
//     double c = (1.0/N) * ((gram_Psi_k(p, p) + gram_l(p, p))*sum_trace/(d/m_k) - d/m_k);
//     Psi_k(p, p) = (-b + sqrt(pow(b, 2.0) - 4.0*a*c))/(2.0*a);
//   }
//   
//   return Psi_k;
// }

// [[Rcpp::export]]
Eigen::MatrixXd update_off_diag(int N, int m_k, double lambda, const Eigen::Map<Eigen::MatrixXd> Psi_k, Eigen::Map<Eigen::MatrixXd> gram_k, Eigen::Map<Eigen::MatrixXd> X_k, Eigen::Map<Eigen::MatrixXd> sum_gram_othermodes, Eigen::Map<Eigen::MatrixXd> W_gram) {
  Eigen::MatrixXd Psi_k_hat = Psi_k;
  
  for (int p = 0; p < m_k - 1; p++){
    for (int q = p + 1; q < m_k; q++){
      Psi_k_hat(p, p) = 0;
      Psi_k_hat(q, q) = 0;
      Psi_k_hat(p, q) = 0;
      Psi_k_hat(q, p) = 0;
      Eigen::MatrixXd gram_Psi_k = gram_k * Psi_k_hat;
      double num = (-1.0/N) * (W_gram(p, q) + W_gram(q, p) + gram_Psi_k(q, p) + gram_Psi_k(p, q) + sum_gram_othermodes(q, p) + sum_gram_othermodes(p, q));
      double denom = (1.0/N) * (gram_k(q, q) + gram_k(p, p));
      Psi_k_hat(p, q) = (num/std::abs(num)) * std::max(std::abs(num) - lambda/N, 0.0)/denom;
      Psi_k_hat(q, p) = Psi_k_hat(p, q);
    }
  }
  
  return Psi_k_hat;
}

// [[Rcpp::export]]
arma::mat stdmvrnormArma(int n, int d) {
  arma::mat Y = arma::randn(n, d);
    return Y;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
Eigen::MatrixXd mvrnormEigen(const Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::SparseMatrix<double> > Omega_sqrt) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::MatrixXd X = solver.compute(Omega_sqrt).solve(Z);
  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXd mat_inv(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd A_inv = A.inverse();
  
  return A_inv;
}

// [[Rcpp::export]]
Eigen::MatrixXd mat_mult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return C;
}

// [[Rcpp::export]]
Eigen::MatrixXd kronecker_prod(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = kroneckerProduct(A, B);

  return C;
}

// [[Rcpp::export]]
Eigen::MatrixXd kronecker_sum(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  int m = A.cols();
  int n = B.cols();
  Eigen::MatrixXd I_m(m, m);
  I_m << Eigen::MatrixXd::Identity(m, m);
  Eigen::MatrixXd I_n(n, n);
  I_n << Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd C = kroneckerProduct(I_n, A) + kroneckerProduct(B, I_m);
  // Eigen::MatrixXd C = kroneckerProduct(A, I_n) + kroneckerProduct(I_m, B);

  return C;
}

