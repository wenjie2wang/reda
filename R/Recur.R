##
## R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
## Copyright (C) 2015-2019
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


##' Recurrent Episodes
##'
##' Specify time segements or recurrent episodes by endpoints.
##'
##' This function is intended to be used for specifying the argument \code{time}
##' in function \code{\link{Recur}}.
##'
##' @name Recur-to
##'
##' @param time1 The left end-points of the recurrent episodes.
##' @param time2 The right end-points of the recurrent episodes.
##'
##' @return A list that consists of two elements named
##'     \code{"time1"} and \code{"time2"}.
##'
##' @export
`%to%` <- function(time1, time2) {
    list(time1 = time1, time2 = time2)
}


##' @rdname Recur-to
##' @export
`%2%` <- `%to%`


##' Formula Response for Recurrent Event Data
##'
##' Create an S4 class object that represents formula response for recurrent
##' event data with optional checking procedures embedded.
##'
##' This is a successor function of the deprecated function \code{Survr}.
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
##' @usage
##'
##' Recur(time, id, event, death, origin, check = TRUE, ...)
##'
##' @param time A numerical vector representing the time of reccurence event or
##'     censoring, or a list with elements named \code{"time1"} and
##'     \code{"time2"} for specifying the follow-up of recurrent events.  In the
##'     latter case, function \code{\%to\%} (or \code{\%2\%}) can be used for
##'     ease of typing.  In addition to numeric values, \code{Date} and
##'     \code{difftime} are allowed and converted to numeric values.  An error
##'     will be thrown if this argument is not specified.
##' @param id Subject identificators. It can be numeric vector, character
##'     vector, or a factor vector.  If it is left unspecified, \code{Recur}
##'     will assume that each row represents a subject.
##' @param event A numeric vector that represents event indicators, event costs,
##'     or event types. Logical vector is allowed and converted to numeric
##'     vector. Censoring is indicated by non-positive values.
##' @param death A numeric vector
##' @param origin The time origin of each subject or process. In addition to
##'     numeric values, \code{Date} and \code{difftime} are also supported and
##'     converted to numeric values.  Different subjects may have different
##'     origins. However, one subject must have the same origin.
##' @param check A logical value suggesting whether to perform data checking
##'     procedure. The default value is \code{TRUE}. \code{FALSE} should be set
##'     with caution and only for processed data already in recerruent event
##'     data framework.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @aliases Recur
##'
##' @return An \code{Recur} object.
##'
##' @export
Recur <- function(time, id, event, death, origin, check = TRUE, ...)
{
    ## warning on `...`
    warn_dots(...)

    ## "time" cannot be missing
    if (missing(time)) {
        stop("The 'time' cannot be missing.")
    }
    time1 <- NULL
    time2 <- process_time(time)
    ## is time1 contained inside of a list of "time"
    if (is.list(time)) {
        time1 <- process_time(time$time1)
        time2 <- process_time(time$time2)
    }
    nRec <- length(time2)

    ## "id" can be left unspecified but cannot contain NA's
    if (missing(id)) {
        ord_id <- sorted_id <- id <- seq_along(time2)
        ID <- factor(id)
        nSubject <- length(time2)
    } else if (anyNA(id)) {
        stop("The 'id' cannot contain missing values.", call. = FALSE)
    } else {
        ## original id's in a character vector
        ID <- factor(id)
        id <- as.numeric(ID)
        uid <- unique(id)
        ord_id <- order(uid)
        sorted_id <- sort(id)
        nSubject <- length(uid)
    }
    first_idx <- ! duplicated(sorted_id)
    last_idx <- ! duplicated(sorted_id, fromLast = TRUE)

    ## "event" can be left unspecified
    ## sort the data by id, time1, time2, and event
    if (! missing(event)) {
        if (isLogicalVector(event)) {
            event <- as.numeric(event)
        }
        ## convert non-positive event all to zero
        ## for consistency with SAS on computing sample MCF
        event[event <= 0] <- 0
        ord <- order(id, time2, - event)
        sorted_event <- event[ord]
    } else {
        ord <- order(id, time2)
        ## all one's before one zero at the end for each id
        sorted_event <- 1 - as.numeric(last_idx)
    }
    sorted_time2 <- time2[ord]

    ## only process origin if 'time1' is not specified
    if (is.null(time1)) {
        ## "origin" is all zero by default
        if (missing(origin)) {
            origin <- 0
        } else {
            if (inherits(origin, "Date")) {
                origin <- as.numeric(origin)
            }
            if (! isNumVector(origin)) {
                stop("The 'origin' must be 'numeric' or 'Date'.",
                     call. = FALSE)
            }
            ## check the length of 'origin'
            len_origin <- length(origin)
            if (len_origin == nRec) {
                origin <- origin[ord][first_idx]
            } else if (len_origin == nSubject) {
                origin <- origin[ord_id]
            } else if (len_origin != 1L) {
                stop(wrapMessages(
                    "The length of 'origin' should be equal to one,",
                    "the number of unique ID's, or the number of records."
                ), call. = FALSE)
            }
        }
        sorted_time1 <- c(NA, sorted_time2[- nRec])
        sorted_time1[which(first_idx)] <- origin
    } else {
        ## throw warning if origin is specified
        if (! missing(origin)) {
            warning(wrapMessages(
                "The specified 'origin' was ignored due to",
                "the presence of 'time1'."
            ), call. = FALSE)
        }
        sorted_time1 <- time1[ord]
    }

    ## "death" can be left unspecified
    ## all censoring by default
    if (missing(death)) {
        sorted_death <- 0
    } else {
        if (isLogicalVector(death)) {
            death <- as.numeric(death)
        }
        ## check the length of 'death'
        len_death <- length(death)
        if (len_death == nRec) {
            sorted_death <- death[ord]
        } else if (len_death == nSubject) {
            sorted_death <- rep(0, nRec)
            sorted_death[which(last_idx)] <- death[ord_id]
        } else if (len_death == 1L) {
            sorted_death <- death
        } else {
            stop(wrapMessages(
                "The length of 'death' should be equal to",
                "the number of unique ID's, or the number of records."
            ), call. = FALSE)
        }
    }

    ## create a numeric matrix for sorted data
    sorted_mat <- cbind(time1 = sorted_time1,
                        time2 = sorted_time2,
                        id = sorted_id,
                        event = sorted_event,
                        death = sorted_death)
    attr(sorted_mat, "ID") <- ID[ord]
    attr(sorted_mat, "ord") <- ord

    ## perform detailed checks
    if (check) {
        check_Recur(sorted_mat, ...)
    }
    ## convert to original ordering
    rev_ord <- order(ord)

    ## return
    methods::new("Recur",
                 .Data = sorted_mat[rev_ord, ],
                 ID = ID,
                 ord = ord,
                 rev_ord = rev_ord,
                 first_idx = which(first_idx),
                 last_idx = which(last_idx),
                 check = check)
}

## helper function that performs detailed checks
check_Recur <- function(sorted_dat, first_idx = NULL, last_idx = NULL)
{
    ## if dat input has an attr 'ID' for internal usage
    sID <- attr(sorted_dat, "ID")
    if (is.null(sID)) {
        sID <- sorted_dat[, "id"]
    }
    sTime1 <- sorted_dat[, "time1"]
    sTime2 <- sorted_dat[, "time2"]
    sEvent <- sorted_dat[, "event"]
    sDeath <- sorted_dat[, "death"]
    sCensor <- sEvent <= 0 | sDeath > 0

    ## set index of the first and the last record of each subject
    if (is.null(first_idx)) {
        first_idx <- ! duplicated(sID)
    }
    if (is.null(last_idx)) {
        last_idx <- ! duplicated(sID, fromLast = TRUE)
    }

    ## stop if event time after censoring time or without censoring time
    idx <- ! sCensor[last_idx]
    if (any(idx)) {
        stop(wrapMessages(
            "Every subject must have one censored time",
            "not earlier than any event time.",
            "Please check subject:",
            paste0(paste(sID[last_idx][idx], collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## stop if more than one censoring time
    cenID <- sID[sCensor]
    idx <- duplicated(cenID)
    if (any(idx)) {
        stop(wrapMessages(
            "Every subject must have only one censored time.",
            "Please check subject:",
            paste0(paste(cenID[idx], collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## stop if missing value of 'time'
    idx <- is.na(sTime1) | is.na(sTime2)
    if (any(idx)) {
        stop(wrapMessages(
            "Event or censoring times cannot be missing.",
            "Please check subject:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## 'time2' has to be later than the 'time1'
    idx <- sTime2 < sTime1
    if (any(idx)) {
        stop(wrapMessages(
            "Event times cannot be earlier than the origin time.",
            "Please check subject:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## return NULL invisibly
    invisible(NULL)
}

## helper function to process 'time'
process_time <- function(x) {
    ## skip if x is null or a list
    if (is.null(x) || is.list(x)) return(x)
    if (inherits(x, "difftime") || inherits(x, "Date"))
        x <- as.numeric(x)
    if (! is.numeric(x))
        stop("The time variables must be 'numeric', 'difftime' or 'Date'.",
             call. = FALSE)
    ## return
    x
}
