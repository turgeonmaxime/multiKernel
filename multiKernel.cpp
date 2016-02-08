#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export("multiKernel_cpp")]]
SEXP multiKernel(SEXP Ys, SEXP Zs, SEXP Ks, double tau) {
  Rcpp::NumericMatrix Yr(Ys);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericMatrix Zr(Zs);                 // creates Rcpp matrix from SEXP
  Rcpp::NumericMatrix Kr(Ks);
  int n = Yr.nrow(), p = Yr.ncol(), q = Zr.ncol();
  
  arma::mat K(Kr.begin(), n, n, false);       // reuses memory and avoids extra copy
  arma::mat Y(Yr.begin(), n, p, false);
  arma::mat Z(Zr.begin(), n, q, false);
  
  // Initialize matrices
  arma::mat Q, R;
  arma::mat B(q, p, fill::ones);
  arma::mat alpha_mat(n, p, fill::zeros);
  arma::vec tuning_vect(n);
  tuning_vect.fill(tau);
  arma::mat tuning_mat(n, n, fill::zeros);
  tuning_mat.diag() = tuning_vect;
  arma::mat weight_mat = inv(K + tuning_mat);
  
  for(int i = 0; i < 1000; i++) {
    // Update B using QR-decomposition
    qr_econ(Q, R, Z);
    B = inv(R) * trans(Q) * (Y - K * alpha_mat);

    // Update alpha_mat
    for(int j = 0; j < p; j++) {
      alpha_mat.col(j) = weight_mat * (Y.col(j) - Z * B.col(j));
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_mat, 
    Rcpp::Named("B") = B) ;
}