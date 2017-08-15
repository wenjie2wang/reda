################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015-2017
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
    origin <- dat[, "origin"]

    if (any(is.na(ID)))
        stop("'ID' cannot contain missing values.")
    if (inherits(time, "difftime") || inherits(time, "Date"))
        time <- unclass(time)
    if (! is.numeric(time))
        stop("Time variable must be 'numeric', 'difftime' or 'Date'.")
    if (any(! event %in% c(0, 1)))
        stop("'event' must be coded as 0 (censoring) or 1 (event).")
    if (inherits(origin, "Date"))
        origin <- unclass(origin)
    if (! is.numeric(origin))
        stop("Origin variable must be 'numeric' or 'Date'.")

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
                   as.data.frame(dat[(ord <- order(ID, time, 1 - event)), ])
               } else {
                   as.data.frame(dat[ord, ])
               }

    if (check) {
        sID <- sortDat[, "ID"]
        sTime <- sortDat[, "time"]
        sEvent <- sortDat[, "event"]
        sOrigin <- sortDat[, "origin"]
        sIDnam <- IDnam[ord]

        ## issue 1: event time after censoring time or without censoring time
        idx <- ! duplicated(sID, fromLast = TRUE) & sEvent != 0
        if (any(idx)) {
            stop("Every subject must have one censored time ",
                 "not earlier than any event time.",
                 "\nPlease check subject: ",
                 paste(sIDnam[idx], collapse = ", "), ".")
        }

        ## issue 2: more than one censoring time
        cenID <- sIDnam[sEvent != 1]
        idx <- duplicated(cenID)
        if (any(idx)) {
            stop("Every subject must have only one censored time.",
                 "\nPlease check subject: ",
                 paste(cenID[idx], collapse = ", "), ".")
        }

        ## stop if missing value of 'time'
        idx <- is.na(sTime)
        if (any(idx)) {
            tmpID <- unique(sIDnam[idx])
            stop("Event or censoring times cannot be missing.",
                 "\nPlease check subject: ",
                 paste(tmpID, collapse = ", "), ".")
        }

        ## 'time' has to be later than the 'origin'
        idx <- sTime < sOrigin
        if (any(idx)) {
            tmpID <- unique(sIDnam[idx])
            stop("Event times cannot be earlier than the origin time.",
                 "\nPlease check subject: ",
                 paste(tmpID, collapse = ", "), ".")
        }

        ## For one subject, the 'origin' has to be the same
        tmp <- tapply(sOrigin, sIDnam, function(a) {
            length(unique(a))
        })
        idx <- tmp > 1
        if (any(idx)) {
            tmpID <- unique(names(tmp)[idx])
            stop("The origin variable has to be the same for one subject.",
                 "\nPlease check subject: ",
                 paste(tmpID, collapse = ", "), ".")
        }
    }
    ## return
    mat <- as.matrix(dat)
    attr(mat, "ID_") <- IDnam
    attr(mat, "check") <- check
    attr(mat, "ord") <- ord
    invisible(mat)
}
