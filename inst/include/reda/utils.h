#ifndef UTILS_H
#define UTILS_H

#include <algorithm>            // std::max, std::set_union, etc.
#include <cmath>                // std::pow and std::sqrt, etc.
#include <limits>
#include <map>
#include <math.h>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <RcppArmadillo.h>


namespace Reda {

    // compare double-precision numbers for almost equality
    inline bool isAlmostEqual(double A, double B)
    {
        double MaxRelDiff {std::numeric_limits<double>::epsilon()};
        // compute the difference.
        double diff = std::abs(A - B);
        A = std::abs(A);
        B = std::abs(B);
        // Find the largest
        double largest = (B > A) ? B : A;
        if (diff <= largest * MaxRelDiff) {
            return true;
        } else {
            return false;
        }
    }
    inline bool is_gt(double A, double B)
    {
        return (! isAlmostEqual(A, B)) && (A > B);
    }
    inline bool is_lt(double A, double B)
    {
        return (! isAlmostEqual(A, B)) && (A < B);
    }
    inline bool is_ge(double A, double B)
    {
        return ! is_lt(A, B);
    }
    inline bool is_le(double A, double B)
    {
        return ! is_gt(A, B);
    }

    // cumulative sum in possibly reverse order
    inline arma::vec cum_sum(const arma::vec& x,
                             const bool reversely = false)
    {
        // if cumsum reversely
        if (reversely) {
            const unsigned long int n_x {x.n_rows};
            arma::vec res {arma::zeros(n_x)};
            double tmp {0.0};
            for (size_t i {1}; i <= n_x; ++i) {
                tmp += x[n_x - i];
                res[n_x - i] = tmp;
            }
            return res;
        }
        // otherwise, using arma::cumsum
        return arma::cumsum(x);
    }

    // aggregate sum of a vector based on same indices
    inline arma::vec aggregate_sum(const arma::vec& x,
                                   const arma::vec& indices,
                                   const bool simplify = true,
                                   const bool cumulative = false,
                                   const bool reversely = false)
    {
        const unsigned int n_x { x.n_elem };
        if (n_x != indices.n_elem) {
            throw std::logic_error(
                "The x and indices must have the same length."
                );
        }
        arma::vec uniInd { arma::unique(indices) };
        const unsigned int n_uniInd { uniInd.n_elem };

        // get sorted x and indices
        arma::uvec sort_ind { arma::sort_index(indices) };
        arma::vec sorted_x = x.elem(sort_ind);
        arma::vec sorted_indices = indices.elem(sort_ind);

        // the x's having a same index are summed
        arma::vec sumVec { arma::zeros(n_uniInd) };

        // early exit if no need to aggregate for all unique indices
        bool is_all_unique { n_uniInd == n_x };
        if (is_all_unique) {
            sumVec = sorted_x;
        } else {
            size_t i {0};
            for (size_t j {0}; j < n_x; ++j) {
                if (! isAlmostEqual(uniInd(i), sorted_indices(j))) {
                    ++i;
                }
                sumVec(i) += sorted_x(j);
            }
        }
        if (cumulative) {
            sumVec = cum_sum(sumVec, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify || is_all_unique) {
            return sumVec;
        }
        // else
        arma::vec out {arma::zeros(n_x)};
        for (size_t i {0}; i < n_x; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {
                    out[i] = sumVec[j];
                    break;
                }
            }
        }
        return out;
    }

    // inline handy functions
    inline arma::vec mat2vec(const arma::mat& x) {
        return arma::conv_to<arma::vec>::from(x);
    }

    // function that returns the indices of the first unique indices
    inline arma::uvec find_first_unique(const arma::uvec& x)
    {
        std::unordered_set<int> seen;
        std::vector<unsigned int> res;
        for (size_t i {0}; i < x.n_rows; ++i) {
                if (seen.insert(x(i)).second) {
                    // if not duplicated, add index to vector res
                    res.push_back(i);
                }
        }
        return arma::conv_to<arma::uvec>::from(res);
    }
    // function that returns the indices of the last unique indices
    inline arma::uvec find_last_unique(const arma::uvec& x)
    {
        std::unordered_set<int> seen;
        std::vector<unsigned int> res;
        for (size_t i {0}; i < x.n_rows; ++i) {
                if (seen.insert(x(x.n_rows - 1 - i)).second) {
                    // if not duplicated, add index to vector res
                    res.push_back(x.n_rows - 1 - i);
                }
        }
        std::reverse(res.begin(), res.end());
        return arma::conv_to<arma::uvec>::from(res);
    }

    // function checking if there exists any duplicates
    inline bool any_duplicated(const arma::vec& x)
    {
        std::unordered_set<double> seen;
        bool res {false};
        for (size_t i {0}; i < x.n_rows; ++i) {
            res = ! seen.insert(x(i)).second;
            if (res) break;
        }
        return res;
    }

    // step function
    inline arma::vec step_fun(const arma::vec& x,
                              const arma::vec& knots,
                              const arma::vec& height)
    {
        // create a map for fast comparison
        std::map<double, double> step_map;
        for (size_t i {0}; i < knots.n_elem; ++i) {
            step_map.insert(std::make_pair(knots(i), height(i + 1)));
        }
        arma::vec res { arma::zeros(x.n_elem) };
        std::map<double, double>::iterator it;
        for (size_t i {0}; i < x.n_elem; ++i) {
            it = step_map.upper_bound(x(i));
            if (it != step_map.begin()) {
                --it;
                res(i) = it->second;
            } else {
                res(i) = height(0);
            }
        }
        return res;
    }

    // quantile function
    // type 5 in quantile
    // reference: Hyndman and Fan (1996)
    inline double arma_quantile(const arma::vec& x, const double prob) {
        const double alpha { 0.5 };
        const unsigned int n { x.n_elem };
        if (prob < (1 - alpha) / n) {
            return x.min();
        }
        if (prob > (n - alpha) / n) {
            return x.max();
        }
        arma::vec inc_x { arma::sort(x) };
        int k { static_cast<int>(std::floor(n * prob + alpha)) };
        double pk { (k - alpha) / n };
        double w { (prob - pk) * n };
        return (1 - w) * inc_x(k - 1) + w * inc_x(k);
    }
    inline arma::vec arma_quantile(const arma::vec& x, const arma::vec& probs) {
        arma::vec res { arma::zeros(probs.n_elem) };
        for (size_t i {0}; i < probs.n_elem; ++i) {
            res(i) = arma_quantile(x, probs(i));
        }
        return res;
    }

    // convert arma vector type to Rcpp vector type
    template <typename T>
    inline Rcpp::NumericVector arma2rvec(const T& x) {
        return Rcpp::NumericVector(x.begin(), x.end());
    }
    // convert Rcpp::NumericVector to arma::colvec
    template <typename T>
    inline arma::vec rvec2arma(const T& x) {
        return arma::vec(x.begin(), x.size(), false);
    }

}

#endif
