#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include<Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector rrisk_exponential(arma::mat z, arma::mat zCoef)
{
  unsigned int nrow_z {z.n_rows};
  Rcpp::NumericVector res(nrow_z);
  for (size_t i {0}; i < nrow_z; ++i) {
    res[i] = std::exp(arma::sum(z.row(i) % zCoef.row(i)));
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector rrisk_linear(arma::mat z, arma::mat zCoef)
{
  unsigned int nrow_z {z.n_rows};
  Rcpp::NumericVector res(nrow_z);
  for (size_t i {0}; i < nrow_z; ++i) {
    res[i] = 1 + arma::sum(z.row(i) % zCoef.row(i));
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector rrisk_excess(arma::mat z, arma::mat zCoef)
{
  unsigned int nrow_z {z.n_rows};
  Rcpp::NumericVector res(nrow_z);
  for (size_t i {0}; i < nrow_z; ++i) {
    res[i] = std::exp(arma::sum(arma::log(1 + z.row(i) % zCoef.row(i))));
  }
  return res;
}
