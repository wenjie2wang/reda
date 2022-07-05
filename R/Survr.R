##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2022
##
## This file is part of the R package reda.
##
## The R package reda is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package reda is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##


## collation after class.R
##' @include class.R
NULL


##' Formula Response for Recurrent Event Data
##'
##' Create an S4 class that represents formula response for recurrent event data
##' modeled by methods based on counts and rate function.  Note that the
##' function is deprecated since version 0.5.0 and will be removed in future.
##'
##' This is a similar function to \code{Survr} in package
##' \pkg{survrec} but with a more considerate checking procedure embedded for
##' recurrent event data modeled by methods based on counts and rate function.
##' The checking rules apply to each subject respectively and include that
##' \itemize{
##'     \item Subject identification, event times, censoring time, and event
##'         indicator cannot be missing or contain missing values.
##'     \item There has to be only one censoring time not earlier than
##'         any event time.
##'     \item The time origin has to be the same and not later than any event
##'         time.
##' }
##'
##' @param ID Subject identificators. It can be numeric vector, character
##'     vector, or a factor vector.
##' @param time Time of reccurence event or censoring. In addition to numeric
##'     values, \code{Date} and \code{difftime} are supported and converted to
##'     numeric values.
##' @param event A numeric vector indicating failure cost or event indicator
##'     taking positive values as costs (\code{1} as events), and non-positive
##'     values as censoring.  Logical vector is allowed and will be converted to
##'     numeric vector.
##' @param origin The time origin of each subject or process. In addition to
##'     numeric values, \code{Date} and \code{difftime} are also supported and
##'     converted to numeric values.  Different subjects may have different
##'     origins. However, one subject must have the same origin.
##' @param check A logical value suggesting whether to perform data checking
##'     procedure. The default value is \code{TRUE}. \code{FALSE} should be set
##'     with caution and only for processed data already in recerruent event
##'     data framework.
##' @param ... Other arguments for future usage.
##'
##' @aliases Survr
##'
##' @export
Survr <- function(ID, time, event, origin = 0, check = TRUE, ...)
{
    ## deprecated from version 0.5.0
    .Deprecated(new = "Recur",
                msg = sprintf(paste(
                    "'%s()' is deprecated and will be removed in future.",
                    "Please use '%s()' instead.",
                    "\n'help(\"Recur\")' for details."),
                    "Survr", "Recur"
                    ))

    if (missing(ID))
        stop("ID variable cannot be missing.")
    if (missing(time))
        stop("Time varibale cannot be missing.")
    if (missing(event))
        stop("Event variable cannot be missing.")

    dat <- data.frame(ID = ID, time = time, event = event, origin = origin)
    check_Survr(dat, check = check)
}


### internal function ==========================================================
check_Survr <- function(dat, check, ...)
{
    ID <- dat[, "ID"]
    time <- dat[, "time"]
    event <- dat[, "event"]
    origin <- dat[, "origin"]

    if (anyNA(ID))
        stop("'ID' cannot contain missing values.", call. = FALSE)
    if (inherits(time, "difftime") || inherits(time, "Date"))
        time <- unclass(time)
    if (! is.numeric(time))
        stop("Time variable must be 'numeric', 'difftime' or 'Date'.",
             call. = FALSE)
    if (inherits(origin, "Date"))
        origin <- unclass(origin)
    if (! isNumVector(origin))
        stop("Origin variable must be 'numeric' or 'Date'.", call. = FALSE)
    if (isLogicalVector(event))
        event <- as.numeric(event)
    ## convert non-positive event all to zero
    ## for ease of later computing sample MCF
    event[event <= 0] <- 0

    ## if dat input has an attr 'ID_' for internal usage
    ID_ <- attr(dat, "ID")
    if (is.null(ID_)) {
        ## check whether 'ID' is numeric or not. convert if not.
        IDnam <- factor(ID)
        dat[, "ID"] <- ID <- as.numeric(IDnam)
    } else {
        IDnam <- ID_
    }

    ## whether data were ordered or not
    ord <- attr(dat, "ord")
    sortDat <- if (is.null(ord) || length(ord) != length(ID_)) {
                   ## sort the data by ID, time, and event
                   as.data.frame(dat[(ord <- order(ID, time, - event)), ])
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
        idx <- ! duplicated(sID, fromLast = TRUE) & sEvent > 0
        if (any(idx)) {
            stop(wrapMessages(
                "Every subject must have one censored time",
                "not earlier than any event time.",
                "Please check subject:",
                paste0(paste(sIDnam[idx], collapse = ", "), ".")
            ), call. = FALSE)
        }

        ## issue 2: more than one censoring time
        cenID <- sIDnam[sEvent <= 0]
        idx <- duplicated(cenID)
        if (any(idx)) {
            stop(wrapMessages(
                "Every subject must have only one censored time.",
                "Please check subject:",
                paste0(paste(cenID[idx], collapse = ", "), ".")
            ), call. = FALSE)
        }

        ## stop if missing value of 'time'
        idx <- is.na(sTime)
        if (any(idx)) {
            tmpID <- unique(sIDnam[idx])
            stop(wrapMessages(
                "Event or censoring times cannot be missing.",
                "Please check subject:",
                paste0(paste(tmpID, collapse = ", "), ".")
            ), call. = FALSE)
        }

        ## stop if missing value of 'origin'
        idx <- is.na(sOrigin)
        if (any(idx)) {
            tmpID <- unique(sIDnam[idx])
            stop(wrapMessages(
                "The origin times cannot be missing.",
                "Please check subject:",
                paste0(paste(tmpID, collapse = ", "), ".")
            ), call. = FALSE)
        }

        ## 'time' has to be later than the 'origin'
        idx <- sTime < sOrigin
        if (any(idx)) {
            tmpID <- unique(sIDnam[idx])
            stop(wrapMessages(
                "Event times cannot be earlier than the origin time.",
                "Please check subject:",
                paste0(paste(tmpID, collapse = ", "), ".")
            ), call. = FALSE)
        }

        ## For one subject, the 'origin' has to be the same
        tmp <- tapply(sOrigin, factor(sIDnam), function(a) {
            length(unique(a))
        })
        idx <- tmp > 1
        if (any(idx)) {
            tmpID <- unique(names(tmp)[idx])
            stop(wrapMessages(
                "The origin variable has to be the same for one subject.",
                "Please check subject:",
                paste0(paste(tmpID, collapse = ", "), ".")
            ), call. = FALSE)
        }
    }
    ## return
    mat <- as.matrix(cbind(ID = ID, time = time,
                           event = event, origin = origin))
    invisible(
        methods::new("Survr",
                     .Data = mat,
                     ID = IDnam,
                     check = check,
                     ord = ord)
    )
}
