#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double computeLeastSq(mat Y, mat K, mat alpha, mat Z, mat B);

// [[Rcpp::export("multiKernel_cpp")]]
SEXP multiKernel(NumericMatrix Yr, NumericMatrix Zr, NumericMatrix Kr, double tau) {
  int n = Yr.nrow(), p = Yr.ncol(), q = Zr.ncol();
  
  arma::mat K(Kr.begin(), n, n, false);       // reuses memory and avoids extra copy
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat Z(Zr.begin(), n, q, false);
  
  // Initialize matrices
  arma::vec tuning_vect(n);
  tuning_vect.fill(tau);
  arma::mat tuning_mat(n, n, fill::zeros);
  tuning_mat.diag() = tuning_vect;
  arma::mat weight_mat = inv(K + tuning_mat);
  
  arma::mat B = inv(trans(Z) * weight_mat * Z) * trans(Z) * weight_mat * Y;
  arma::mat alpha_mat = weight_mat * (Y - Z * B);
  
  double newLS = computeLeastSq(Y, K, alpha_mat, Z, B);
  double BIC = 2 * newLS + log(n) * as_scalar(accu(B > 0) + accu(alpha_mat > 0));
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_mat, 
    Rcpp::Named("B") = B,
    Rcpp::Named("LS") = newLS,
    Rcpp::Named("BIC") = BIC);
}

// [[Rcpp::export("multiKernel_noCon_cpp")]]
SEXP multiKernel_noCon(NumericMatrix Yr, NumericMatrix Kr, double tau) {
  int n = Yr.nrow(), p = Yr.ncol();
  
  arma::mat K(Kr.begin(), n, n, false);       // reuses memory and avoids extra copy
  arma::mat Y(Yr.begin(), n, p, false);

  // Initialize matrices
  arma::vec tuning_vect(n);
  tuning_vect.fill(tau);
  arma::mat tuning_mat(n, n, fill::zeros);
  tuning_mat.diag() = tuning_vect;

  arma::mat alpha_mat = inv(K + tuning_mat) * Y;
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_mat);
}

double computeLeastSq(mat Y, mat K, mat alpha, mat Z, mat B) {
  double res = 0.5 * norm(Y - K * alpha - Z * B, "fro");
  return res;
}