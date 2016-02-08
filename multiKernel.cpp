#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double computeLeastSq(mat Y, mat K, mat alpha, mat Z, mat B);

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
  
  // Initialize loss value
  double currentLS = 1.0;
  double newLS = 0.0;
  int counter = 0;
  
  while(fabs(currentLS - newLS) > 0.00000001 && counter < 10000) {
    counter++;
    currentLS = newLS;
    // Update B using QR-decomposition
    qr_econ(Q, R, Z);
    B = inv(R) * trans(Q) * (Y - K * alpha_mat);

    // Update alpha_mat
    for(int j = 0; j < p; j++) {
      alpha_mat.col(j) = weight_mat * (Y.col(j) - Z * B.col(j));
    }
    
    newLS = computeLeastSq(Y, K, alpha_mat, Z, B);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha_mat, 
    Rcpp::Named("B") = B,
    Rcpp::Named("iter") = counter);
}

double computeLeastSq(mat Y, mat K, mat alpha, mat Z, mat B) {
  double res = - 0.5 * norm(Y - K * alpha - Z * B, "fro");
  
  return res;
}