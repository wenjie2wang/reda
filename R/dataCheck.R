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

    ID <- dat[, "ID"]
    time <- dat[, "time"]
    event <- dat[, "event"]

    ## event time must be numeric
    if (! is.numeric(time))
        stop("'time' must be numeric.")
    ## check missing value on 'ID'
    if (any(is.na(ID)))
        stop("'ID' cannot be missing.")
    ## check coding and missing value on 'event'
    if (any(! event %in% c(0, 1)))
        stop("'event' must be coded as 0 (censoring) or 1 (event).")

    ## if dat input has an attr 'ID_' for internal usage
    ID_ <- attr(dat, "ID_")
    if (is.null(ID_)) {
        ## check whether 'ID' is numeric or not. convert if not.
        IDnam <- factor(ID)
        dat[, "ID"] <- ID <- as.numeric(IDnam)
    } else {
        IDnam <- ID_
    }

    ## whether data were ordered or not
    ord <- attr(dat, "ord")
    sortDat <- if (is.null(ord)) {
                   ## sort the data by ID, time, and event
                   as.data.frame(dat[(ord <- order(ID, time, event)), ])
               } else {
                   as.data.frame(dat[ord, ])
               }

    if (check) {
        sID <- sortDat[, "ID"]
        sTime <- sortDat[, "time"]
        sEvent <- sortDat[, "event"]
        sIDnam <- IDnam[ord]

        ## issue 1: event time after censoring time or without censoring time
        idx1 <- ! duplicated(sID, fromLast = TRUE) & sEvent != 0
        if (any(idx1)) {
            stop(paste("Every subject must have one censored time",
                       "later than event times.",
                       "\nPlease check subject:",
                       paste(sIDnam[idx1], collapse = ", ")))
        }

        ## issue 2: more than one censoring time
        cenID <- sIDnam[sEvent != 1]
        idx2 <- duplicated(cenID)
        if (any(idx2)) {
            stop(paste("Every subject must have only one censored time.",
                       "\nPlease check subject:",
                       paste(cenID[idx2], collapse = ", ")))
        }

        ## stop if missing value of 'time'
        idx3 <- is.na(sTime)
        if (any(idx3)) {
            tmpID <- unique(sIDnam[idx3])
            stop(paste("Event or censoring times cannot be missing.",
                       "\nPlease check subject:",
                       paste(tmpID, collapse = ", ")))
        }
    }
    ## return
    mat <- as.matrix(dat)
    attr(mat, "ID_") <- IDnam
    attr(mat, "check") <- check
    attr(mat, "ord") <- ord
    invisible(mat)
}
