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


##' Formula Response for Recurrent Event Data
##'
##' \code{Survr} is an S3 class that represents
##' formula response for recurrent event data
##' modeled by methods based on counts and rate function.
##' The last letter 'r' in 'Survr' represents 'rate'.
##'
##' This is a similar function to \code{Survr} in package
##' \pkg{survrec} but with a more considerate checking procedure embedded for
##' recurrent event data modeled by methods based on counts and rate function.
##' The checking rules apply to each subject and include that
##' \itemize{
##'     \item Subject identification, event times, censoring time, and event
##'         indicator cannot be missing or contain missing values.
##'     \item Event indicator must be coded as 0 (censored) or 1 (event).
##'     \item There has to be only one censoring time not earlier than
##'         any event time.
##'     \item The time origin has to be the same and not later than any event
##'         time.
##' }
##'
##' @param ID Subject identificator.
##' @param time Time of reccurence event or censoring. In addition to numeric
##'     values, \code{Date} and \code{difftime} are also supported.
##' @param event The event indicator, \code{0 = censored}, \code{1 = event}.
##' @param origin The time origin of each subject or process. Different subjects
##'     may have different origins. However, one subject must have the same
##'     origin.
##' @param check Logical value suggesting whether to perform data checking
##'     procedure. The default value is \code{TRUE}. \code{FALSE} should be set
##'     with caution and only for processed data already in recerruent event
##'     data framework.
##' @param ... Other arguments for future usage.
##' @aliases Survr
##' @seealso \code{\link{rateReg}} for model fitting.
##' @export
Survr <- function(ID, time, event, origin = 0, check = TRUE, ...) {
    if (missing(ID))
        stop("ID variable cannot be missing.")
    if (missing(time))
        stop("Time varibale cannot be missing.")
    if (missing(event))
        stop("Event variable cannot be missing.")

    dat <- data.frame(ID = ID, time = time, event = event, origin = origin)
    dat <- check_Survr(dat, check = check)
    attr(dat, "check") <- check
    class(dat) <- c("matrix", "Survr")
    invisible(dat)
}


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
