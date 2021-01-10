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

#include <Rcpp.h>

// modified from https://gallery.rcpp.org/articles/fast-factor-generation/
// cannot handle NA's
template <int RTYPE>
Rcpp::List rcpp_factorize_template(const Rcpp::Vector<RTYPE>& x)
{
    Rcpp::Vector<RTYPE> levs { Rcpp::sort_unique(x) };
    Rcpp::IntegerVector out { Rcpp::match(x, levs) };
    return Rcpp::List::create(
        Rcpp::Named("ID") = Rcpp::as<Rcpp::CharacterVector>(levs),
        Rcpp::Named("id") = out
        );
}

// [[Rcpp::export]]
Rcpp::List rcpp_factorize(SEXP x)
{
    switch(TYPEOF(x)) {
        case INTSXP:
            return rcpp_factorize_template<INTSXP>(x);
        case REALSXP:
            return rcpp_factorize_template<REALSXP>(x);
        case STRSXP:
            return rcpp_factorize_template<STRSXP>(x);
    }
    return Rcpp::List();
}
