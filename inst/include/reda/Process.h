#ifndef PROCESS_H
#define PROCESS_H

#include <RcppArmadillo.h>
#include "utils.h"

namespace Reda {

    // Process class
    class Process {
    private:
        unsigned int id;
        arma::vec time1;
        arma::vec time2;
        arma::vec event;
        double origin_time;
        double censor_time;

    public:
        // constructors
        Process();

        Process(unsigned int id_, arma::vec time1_,
                arma::vec time2_, arma::vec event_) :
            id(id_), time1(time1_), time2(time2_), event(event_)
        {
            censor_time = arma::max(time2);
            origin_time = arma::min(time1);
            // sort by time1, time2 and - event
            arma::uvec sorted_idx { arma::sort_index(time1 + time2 - event) };
            time1 = time1.elem(sorted_idx);
            time2 = time2.elem(sorted_idx);
            event = event.elem(sorted_idx);
        }

        // get function members
        inline unsigned int get_id() { return id; }
        inline arma::vec get_time1() { return time1; }
        inline arma::vec get_time2() { return time2; }
        inline arma::vec get_event() { return event; }
        inline double get_origin_time() { return origin_time; }
        inline double get_censor_time() { return censor_time; }

        // set function members
        inline void set_id(unsigned int id_) { id = id_; }

        // other method members
        inline int is_at_risk(double time)
        {
            if (is_gt(time, censor_time) ||
                is_gt(origin_time, time)) {
                // early exit
                return 0;
            }
            for (size_t i {0}; i < time1.n_rows; ++i) {
                if (is_ge(time, time1(i)) && is_le(time, time2(i))) {
                    return 1;
                }
            }
            return 0;
        }


    };                          // end of Process class


}


#endif
