//
// R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
// Copyright (C) 2015-2021
//
// This file is part of the R package reda.
//
// The R package reda is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package reda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


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
