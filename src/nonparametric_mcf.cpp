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

#include <RcppArmadillo.h>
#include <reda.h>

// [[Rcpp::export]]
Rcpp::List cpp_np_mcf(
    const arma::vec& time1,
    const arma::vec& time2,
    const arma::uvec& id,
    const arma::vec& event,
    const unsigned int& point_method = 1,
    const unsigned int& var_method = 1,
    const unsigned int& ci_method = 1,
    const double& ci_level = 0.95,
    const unsigned int& var_bootstrap_method = 1,
    const unsigned int& var_bootstrap_B = 30
    )
{
    Reda::MCF mcf_obj { Reda::MCF(time1, time2, id, event) };
    mcf_obj.estimate(
        point_method,
        var_method,
        ci_method,
        ci_level,
        var_bootstrap_method,
        var_bootstrap_B
        );
    return Rcpp::List::create(
        Rcpp::Named("jump_time") = mcf_obj.jump_time,
        Rcpp::Named("riskset_size") = mcf_obj.riskset_size,
        Rcpp::Named("inst_rate") = mcf_obj.inst_rate,
        Rcpp::Named("cum_rate") = mcf_obj.cum_rate,
        Rcpp::Named("se_cum_rate") = mcf_obj.se_cum_rate,
        Rcpp::Named("lower_cum_rate") = mcf_obj.lower_cum_rate,
        Rcpp::Named("upper_cum_rate") = mcf_obj.upper_cum_rate
        );
}
