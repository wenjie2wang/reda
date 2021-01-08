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

#ifndef NONPARAMETRIC_MCF_H
#define NONPARAMETRIC_MCF_H

#include <vector>
#include <RcppArmadillo.h>
#include "Process.h"
#include "utils.h"

namespace Reda {

    class MCF {
    private:
        arma::uvec s_id_;       // sorted process ID's
        arma::uvec u_id_;       // sorted unique process ID's
        arma::vec s_time1_;     // time 1 of recurrent episode sorted by ID
        arma::vec s_time2_;     // time 2 of recurrent episode sorted by ID
        arma::vec s_event_;     // events sorted by ID
        // a vector of processes
        std::vector<Process> processes_;

        // 0: no estimates by default
        // 1: for Nelsen-Aalen estimator
        // 2: for Cook-Lawless estimator (csm)
        unsigned int point_method_ {0};

        // method for variance estimates of MCF
        // 0: no estimates by default
        // 1: Lawless and Nadeau method (1995)
        // 2: Poisson process method
        // 3: bootstrap method
        // 4: Cook and Lawless's cumulative sample variance (2007)
        unsigned int var_method_ {0};

        // 0: no ci estimates
        // 1: normality of mcf estimates
        // 2: normality of log(mcf estimates)
        // 3: quantile of mcf estimates from bootstrap samples
        unsigned int ci_method_ {0};
        double ci_level_ {0.95};

        // method for variance estimates from bootstrap samples
        // 0: no variance estimates
        // 1: sd of mcf estimates from bootstrap samples
        // 2: quantile of mcf estimates based on normality
        unsigned int var_bootstrap_method_ {0};
        // number of bootstrap samples
        unsigned int var_bootstrap_B_ {0};
        // upper quantile
        double upper_quantile_ {0.75};
        double lower_quantile_ {0.25};

        // private methods without much care
        void compute_point_estimate();
        void compute_var_ln();
        void compute_var_poisson();
        arma::mat compute_var_bootstrap();
        void compute_var_csv();

    public:
        double min_origin;       // the earliest origin time
        arma::vec jump_time;     // unique jump times where inst_rate > 0

        arma::vec inst_rate;     // instant hazard rate at jump_time
        arma::vec cum_rate;      // cumulative hazard rate at jump_time
        arma::uvec riskset_size; // risk-set size at t-0

        // variance estimates of cumulative hazard rate at jump_time
        arma::vec var_cum_rate;
        arma::vec se_cum_rate;

        // confidence interval estimates
        arma::vec lower_cum_rate;
        arma::vec upper_cum_rate;

        // default constructor
        MCF() {}

        // constructors
        MCF(const arma::vec& time1,
            const arma::vec& time2,
            const arma::uvec& id,
            const arma::vec& event)
        {
            arma::uvec id_sort_idx { arma::sort_index(id, "ascend") };
            s_id_ = id.elem(id_sort_idx);
            s_time1_ = time1.elem(id_sort_idx);
            s_time2_ = time2.elem(id_sort_idx);
            s_event_ = event.elem(id_sort_idx);
            u_id_ = arma::unique(id); // sorted unique id
            arma::uvec id_first_idx { find_first_unique(s_id_) };
            arma::uvec id_last_idx { find_last_unique(s_id_) };
            // create processes
            for (size_t i {0}; i < u_id_.n_rows; ++i) {
                arma::uvec idx {
                    arma::regspace<arma::uvec>(id_first_idx(i), id_last_idx(i))
                };
                processes_.push_back(
                    Process(u_id_(i),
                            s_time1_.elem(idx),
                            s_time2_.elem(idx),
                            s_event_.elem(idx))
                    );
            }
        }

        // function members
        // main user interfance
        void estimate(
            const unsigned int& point_method,
            const unsigned int& var_method,
            const unsigned int& ci_method,
            const double& ci_level,
            // ignored unless var_method = 3
            const unsigned int& var_bootstrap_method,
            const unsigned int& var_bootstrap_B
            );

    };                          // end of class definition of MCF


    // point estimates of MCF
    inline void MCF::compute_point_estimate()
    {
        // locate unique jump times
        jump_time = arma::unique(s_time2_);
        riskset_size = arma::zeros<arma::uvec>(jump_time.n_rows);
        // compute risk-set size based on point_method
        switch(point_method_) {
            default:
            case 1: {
                // get the size of risk set at jump times
                for (size_t i {0}; i < jump_time.n_rows; ++i) {
                    for (size_t j {0}; j < processes_.size(); ++j) {
                        riskset_size(i) +=
                            processes_.at(j).is_at_risk(jump_time(i));
                    }
                }
                break;
            }
            case 2: {
                // the riskset size is a constant for CSM
                riskset_size = riskset_size.fill(u_id_.n_rows);
                break;
            }
        }
        // get the jump size
        arma::vec jump_size { aggregate_sum(s_event_, s_time2_) };
        // get the minimum time origin
        min_origin = arma::min(s_time1_);
        // compute the instant hazard_ rate
        inst_rate = jump_size / riskset_size;
        cum_rate = cum_sum(inst_rate);
    }
    // compute variance estimates by Lawless and Nadeau (1995)
    inline void MCF::compute_var_ln()
    {
        // create aliases for ease of following notations in the paper
        const arma::vec& tj = jump_time;
        const arma::uvec& delta_tj = riskset_size;
        const arma::vec& m_tj = inst_rate;

        // initializations
        var_cum_rate = arma::zeros(jump_time.n_rows);
        // for each process
        for (size_t i {0}; i < processes_.size(); ++i) {
            arma::vec res_i { arma::zeros(jump_time.n_rows) };
            arma::vec tj_i { processes_.at(i).get_time2() };
            arma::vec event_i { processes_.at(i).get_event() };
            // expand results to jump_time
            tj_i = arma::join_cols(tj_i, tj);
            event_i = arma::join_cols(event_i,
                                      arma::zeros(jump_time.n_rows));
            arma::vec n_i_tj {
                aggregate_sum(event_i, tj_i)
            };
            // for each jump time
            for (size_t j {0}; j < jump_time.n_rows; ++j) {
                if (isAlmostEqual(delta_tj(j), 0.0)) {
                    res_i(j) = 0;
                }
                // compute at-risk indicator
                double delta_i_tj { static_cast<double>(
                        processes_.at(i).is_at_risk(jump_time(j))
                        ) };
                res_i(j) = delta_i_tj / delta_tj(j) * (n_i_tj(j) - m_tj(j));
            }
            var_cum_rate += arma::pow(cum_sum(res_i), 2);
        }
        se_cum_rate = arma::sqrt(var_cum_rate);
    }


    // compute variance estimates of MCF by Poisson model
    inline void MCF::compute_var_poisson()
    {
        var_cum_rate = arma::zeros(jump_time.n_rows);
        for (size_t i {0}; i < jump_time.n_rows; ++i) {
            if (riskset_size(i) > 0) {
                var_cum_rate(i) = inst_rate(i) / riskset_size(i);
            }
        }
        var_cum_rate = cum_sum(var_cum_rate);
        se_cum_rate = arma::sqrt(var_cum_rate);
    }

    // compute variance estimates of MCF by bootstrap method
    inline arma::mat MCF::compute_var_bootstrap()
    {
        // initialize bootstrap matrix
        arma::mat boot_mat {
            arma::zeros(jump_time.n_rows, var_bootstrap_B_)
                };
        // main loop
        for (size_t b {0}; b < var_bootstrap_B_; ++b) {
            // re-sample processes_ with replacement
            arma::uvec boot_ind;
            if (var_bootstrap_B_ > 0) {
                arma::vec tmp { arma::randu(u_id_.n_rows) };
                tmp = arma::floor(tmp * u_id_.n_rows);
                boot_ind = arma::conv_to<arma::uvec>::from(tmp);
            }
            std::vector<Process> boot_processes { processes_ };
            // get the bootstrap processes_
            arma::vec boot_time2, boot_event;
            for (size_t i {0}; i < u_id_.n_rows; ++i) {
                Process tmp { processes_.at(boot_ind(i)) };
                tmp.set_id(i);
                boot_processes.at(i) = tmp;
                boot_time2 = arma::join_cols(
                    boot_time2, tmp.get_time2()
                    );
                boot_event = arma::join_cols(
                    boot_event, tmp.get_event()
                    );
            }
            // compute point estimates of MCF
            // unique jump times
            arma::vec boot_jump_time { arma::unique(boot_time2) };
            arma::uvec boot_riskset_size {
                arma::zeros<arma::uvec>(boot_jump_time.n_rows)
            };
            // compute risk-set size based on point_method
            switch(point_method_) {
                default:
                case 1: {
                    // get the size of risk set at jump times
                    for (size_t j {0}; j < boot_jump_time.n_rows; ++j) {
                        for (size_t i {0}; i < boot_processes.size(); ++i) {
                            boot_riskset_size(j) +=
                                boot_processes.at(i).is_at_risk(
                                    boot_jump_time(j)
                                    );
                        }
                    }
                    break;
                }
                case 2: {
                    // the riskset size is a constant for CSM
                    boot_riskset_size = boot_riskset_size.fill(u_id_.n_rows);
                    break;
                }
            }
            // get the jump size
            arma::vec boot_jump_size {
                aggregate_sum(boot_event, boot_time2)
            };
            // compute the instant hazard_ rate
            arma::vec boot_inst_rate { boot_jump_size / boot_riskset_size };
            // use step function
            // to make sure bootstrap estimates are at the same time grid
            arma::vec boot_cum_rate {
                arma::join_cols(arma::zeros(1), cum_sum(boot_inst_rate))
            };
            arma::vec boot_res {
                step_fun(jump_time, boot_jump_time, boot_cum_rate)
            };
            boot_mat.col(b) = boot_res;
        } // end of the main loop
        // compute variance and se estimates
        switch(var_bootstrap_method_) {
            case 0:
            default:
                break;
            case 1: {
                var_cum_rate = arma::var(boot_mat, 0, 1);
                se_cum_rate = arma::sqrt(var_cum_rate);
                break;
            }
            case 2: {
                se_cum_rate = arma::zeros(jump_time.n_rows);
                double q_norm_val {
                    R::qnorm(upper_quantile_, 0, 1, 1, 0) -
                        R::qnorm(lower_quantile_, 0, 1, 1, 0)
                        };
                for (size_t j {0}; j < boot_mat.n_rows; ++j) {
                    arma::vec tmp { boot_mat.row(j).t() };
                    double q1 { arma_quantile(tmp, lower_quantile_) };
                    double q2 { arma_quantile(tmp, upper_quantile_) };
                    se_cum_rate(j) = (q2 - q1) / q_norm_val;
                }
                var_cum_rate = arma::pow(se_cum_rate, 2);
                break;
            }
        }
        return boot_mat;
    }

    inline void MCF::compute_var_csv()
    {
        // initializations
        var_cum_rate = arma::zeros(jump_time.n_rows);
        // for each process
        for (size_t i {0}; i < processes_.size(); ++i) {
            arma::vec tj_i { processes_.at(i).get_time2() };
            arma::vec event_i { processes_.at(i).get_event() };
            // expand results to jump_time
            tj_i = arma::join_cols(tj_i, jump_time);
            event_i = arma::join_cols(event_i,
                                      arma::zeros(jump_time.n_rows));
            arma::vec cum_n_i_tj {
                aggregate_sum(event_i, tj_i, true, true)
            };
            arma::vec res_i {
                arma::pow(cum_n_i_tj - cum_rate, 2) / (u_id_.n_rows - 1)
            };
            var_cum_rate += res_i;
        }
        var_cum_rate /= u_id_.n_rows;
        se_cum_rate = arma::sqrt(var_cum_rate);
    }

    inline void MCF::estimate(
        const unsigned int& point_method = 1,
        const unsigned int& var_method = 1,
        const unsigned int& ci_method = 1,
        const double& ci_level = 0.95,
        // ignored unless var_method = 3
        const unsigned int& var_bootstrap_method = 1,
        const unsigned int& var_bootstrap_B = 30
        )
    {
        // point estimates
        this->point_method_ = point_method;
        if (point_method == 0) return;
        this->compute_point_estimate();
        // variance estimates of MCF
        this->var_method_ = var_method;
        arma::mat boot_mat;
        switch(var_method) {
            // if no variance estimates
            case 0:
                return;
            default:
            case 1:
                this->compute_var_ln();
                break;
            case 2:
                this->compute_var_poisson();
                break;
            case 3: {
                this->var_bootstrap_method_ = var_bootstrap_method;
                this->var_bootstrap_B_ = var_bootstrap_B;
                boot_mat = this->compute_var_bootstrap();
                break;
            }
            case 4:
                this->compute_var_csv();
                break;
        }
        // compute confidence interval estimates of MCF
        this->ci_method_ = ci_method;
        this->ci_level_ = ci_level;
        double upper_q { (ci_level + 1) / 2 };
        double cri_val { R::qnorm(upper_q, 0, 1, 1, 0) };
        arma::vec cri_vec { cri_val * se_cum_rate };
        switch(ci_method) {
            case 0:
            default:
                return;
            case 1:
                this->lower_cum_rate = cum_rate - cri_vec;
                this->upper_cum_rate = cum_rate + cri_vec;
                break;
            case 2: {
                arma::vec w_exp { cri_vec };
                for (size_t j {0}; j < cri_vec.n_rows; ++j) {
                    if (is_gt(cum_rate(j), 0)) {
                        w_exp(j) = cri_vec(j) / cum_rate(j);
                    } else {
                        w_exp(j) = 0;
                    }
                }
                w_exp = arma::exp(w_exp);
                this->upper_cum_rate = cum_rate % w_exp;
                this->lower_cum_rate = cum_rate / w_exp;
                break;
            }
            case 3: {
                this->lower_cum_rate = arma::zeros(jump_time.n_rows);
                this->upper_cum_rate = arma::zeros(jump_time.n_rows);
                for (size_t j {0}; j < jump_time.n_rows; ++j) {
                    arma::vec tmp { boot_mat.row(j).t() };
                    this->lower_cum_rate(j) = arma_quantile(tmp, 1 - upper_q);
                    this->upper_cum_rate(j) = arma_quantile(tmp, upper_q);
                }
                break;
            }
        }
    } // end of estimate


}




#endif
