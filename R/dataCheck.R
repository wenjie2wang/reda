################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2016
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
##' @include class.R
NULL


### internal function ==========================================================
check_Survr <- function(dat, check, ...) {

    ## nonsense, just to suppress Note from R CMD check --as-cran
    event <- NULL

    ## check missing value on 'ID'
    if (any(is.na(dat$ID)))
        stop("'ID' cannot be missing.")

    ## check coding and missing value on 'event'
    if (any(! dat$event %in% c(0, 1)))
        stop("'event' must be coded as 0 (censoring) or 1 (event).")

    ## sort the data by ID, time, and event
    dat <- dat[(ord <- with(dat, order(ID, time, event))), ]

    ## if dat input has an attr 'ID_'
    nID <- attr(dat, "ID_")
    if (is.null(nID)) {
        ## check whether 'ID' is numeric or not. convert if not.
        dat$IDnam <- factor(dat$ID, levels = unique(dat$ID))
        dat$ID <- as.numeric(dat$IDnam)
    } else {
        dat$IDnam <- nID
    }

    if (check) {
        ## issue 1: event time after censoring time or without censoring time
        idx1 <- ! duplicated(dat$ID, fromLast = TRUE) & dat$event != 0
        if (any(idx1)) {
            warning(paste("Every subject must have one censored time",
                          "later than event times."))
            stop(paste("Please check subject:",
                       paste(dat$IDnam[idx1], collapse = ", ")))
        }

        ## issue 2: more than one censoring time
        cenID <- subset(dat, event != 1)[, "IDnam"]
        idx2 <- duplicated(cenID)
        if (any(idx2)) {
            warning("Every subject must have only one censored time.")
            stop(paste("Please check subject:",
                       paste(cenID[idx2], collapse = ", ")))
        }

        ## stop if missing value of 'time'
        idx3 <- is.na(dat$time)
        if (any(idx3)) {
            tmpID <- unique(dat$IDnam[idx3])
            warning("Event or censoring times cannot be missing.")
            stop(paste("Please check subject:", paste(tmpID, collapse = ", ")))
        }
    }

    ## return
    attr(dat, "ID_") <- dat$IDnam
    dat$IDnam <- NULL
    mat <- as.matrix(dat)
    attr(mat, "ID_") <- attr(dat, "ID_")
    attr(mat, "ord") <- ord
    attr(mat, "check") <- check
    invisible(mat)
}
