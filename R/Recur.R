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
##' This is a successor function of the deprecated function \code{Survr}.  See
##' the vignette by `vignette("reda-Recur")` for details.
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
##' @param event A numeric vector that may represent the status, costs, or types
##'     of the recurrent events. Logical vector is allowed and converted to
##'     numeric vector. Non-positive values are internally converted to zero
##'     indicating censoring status.
##' @param terminal A numeric vector that may represent the status, costs, or
##'     types of the terminal events.  Logival vector is allowed and converted
##'     to numeric vector.  Non-positive values are internally converted to zero
##'     indicating censoring status.  If a scalar value is specified, all
##'     subjects will have the same status of terminal events at their last
##'     recurrent episodes.  The length of the specified \code{terminal} should
##'     be equal to the number of subjects, or number of data rows.  In the
##'     latter case, each subject may have at most one positive entry of
##'     \code{terminal} at the last recurrent episode.
##' @param origin The time origin of each subject.  If a scalar value is
##'     specified, all subjects will have the same origin at the specified
##'     value.  The length of the specified \code{origin} should be equal to the
##'     number of subjects, or number of data rows.  In the latter case,
##'     different subjects may have different origins.  However, one subject
##'     must have the same origin.  In addition to numeric values, \code{Date}
##'     and \code{difftime} are also supported and converted to numeric values.
##' @param check A character value specifying how to perform the checks for
##'     recurrent event data.  Errors or warnings will be thrown, respectively,
##'     if the \code{check} is specified to be \code{"hard"} (by default) or
##'     \code{"soft"}.  If \code{check = "none"} is specified, no data checking
##'     procedure will be run.
##' @param ... Other arguments for future usage.  A warning will be thrown if
##'     any invalid argument is specified.
##'
##' @aliases Recur
##'
##' @return An \code{Recur} object.
##'
##' @example inst/examples/ex_Recur.R
##'
##' @export
Recur <- function(time, id, event, terminal, origin,
                  check = c("hard", "soft", "none"), ...)
{
    ## warning on `...`
    warn_dots(...)

    ## get check; make it faster by specifying choices
    check <- match.arg(check, choices = c("hard", "soft", "none"))

    ## "time" cannot be missing
    if (missing(time)) {
        stop("The 'time' cannot be missing.")
    }

    ## is time1 contained inside of a list of "time"
    if (is.list(time)) {
        time1 <- process_time(time$time1, allow_null = TRUE)
        time2 <- process_time(time$time2, allow_null = FALSE)
        cl1 <- class(time$time1)
        cl2 <- class(time$time2)
        if (! isTRUE(all.equal.character(cl1, cl2, check.attributes = FALSE))) {
            warning("The 'time1' and 'time2' do not have the same class.")
        }
        time_class <- union(cl1, cl2)
    } else {
        time1 <- NULL
        time2 <- process_time(time, allow_null = FALSE)
        time_class <- class(time)
    }
    nRec <- length(time2)

    ## "id" can be left unspecified but cannot contain NA's
    if (missing(id)) {
        first_idx <- last_idx <- ord_id <- sorted_id <-
            id <- uid <- seq_len(nRec)
        ID <- as.character(id)
        nSubject <- nRec
    } else if (anyNA(id)) {
        stop("The 'id' cannot contain missing values.", call. = FALSE)
    } else {
        ## original id's in a character vector
        id_list <- rcpp_factorize(id)
        ID <- id_list$ID
        id <- id_list$id
        ## ID[id] = input id
        uid <- unique(id)
        ord_id <- order(uid)
        sorted_id <- sort(id)
        nSubject <- length(uid)
        first_idx <- which(! duplicated(sorted_id))
        last_idx <- which(! duplicated(sorted_id, fromLast = TRUE))
    }

    ## "event" can be left unspecified
    ## sort the data by id, time1, time2, and event
    if (! missing(event)) {
        if (isLogicalVector(event, error_na = TRUE)) {
            event <- as.numeric(event)
        } else if (! isNumVector(event, error_na = TRUE)) {
            stop("Invalid 'event', see '?Recur' for details.",
                 call. = FALSE)
        }
        ## convert non-positive event all to zero
        ## for consistency with SAS on computing sample MCF
        event[event <= 0] <- 0
        ## sort by id, time1, (time2,) and - event
        ord <- order(id, time2, - event)
        sorted_event <- event[ord]
    } else {
        ord <- order(id, time2)
        ## all one's before one zero at the end for each id
        sorted_event <- rep.int(1, nRec)
        sorted_event[last_idx] <- 0
    }
    sorted_time2 <- time2[ord]

    ## only process origin if 'time1' is not specified
    if (is.null(time1)) {
        ## "origin" is all zero by default
        if (missing(origin)) {
            sorted_origin <- origin <- 0
        } else {
            origin <- process_time(origin, allow_null = TRUE)
            ## check the length of 'origin'
            len_origin <- length(origin)
            if (len_origin == nRec) {
                sorted_origin <- origin[ord]
                origin <- sorted_origin[first_idx]
            } else if (len_origin == nSubject) {
                origin <- origin[ord_id]
                sorted_origin <- origin[sorted_id]
            } else if (len_origin == 1L) {
                sorted_origin <- origin
            } else {
                stop(wrapMessages(
                    "Invalid length for 'origin'. See '?Recur' for details."
                ), call. = FALSE)
            }
        }
        sorted_time1 <- c(NA, sorted_time2[- nRec])
        sorted_time1[first_idx] <- origin
    } else {
        ## throw warning if origin is specified
        if (! missing(origin)) {
            warning(wrapMessages(
                "The specified 'origin' was ignored due to given 'time1'."
            ), call. = FALSE)
        }
        sorted_time1 <- time1[ord]
        origin <- sorted_time1[first_idx]
        ## sorted_id can be also used as index
        sorted_origin <- origin[sorted_id]
    }

    ## "terminal" can be left unspecified
    if (missing(terminal)) {
        ## all censoring by default
        sorted_terminal <- rep(0, nRec)
    } else {
        if (isLogicalVector(terminal, error_na = TRUE)) {
            terminal <- as.numeric(terminal)
        } else if (! isNumVector(terminal, error_na = TRUE)) {
            stop("Invalid 'terminal'.  See '?Recur' for details.",
                 call. = FALSE)
        }
        ## check the length of 'terminal'
        len_terminal <- length(terminal)
        if (len_terminal == nRec) {
            sorted_terminal <- terminal[ord]
        } else {
            sorted_terminal <- rep(0, nRec)
            if (len_terminal == nSubject) {
                sorted_terminal[last_idx] <- terminal[ord_id]
            } else if (len_terminal == 1L) {
                sorted_terminal[last_idx] <- terminal
            } else {
                stop(wrapMessages(
                    "Invalid length for 'terminal'.  See '?Recur' for details."
                ), call. = FALSE)
            }
        }
    }

    ## create a numeric matrix for sorted data
    sorted_mat <- cbind(time1 = sorted_time1,
                        time2 = sorted_time2,
                        id = sorted_id,
                        event = sorted_event,
                        terminal = sorted_terminal,
                        origin = sorted_origin)

    ## convert to original ordering
    rev_ord <- order(ord)

    ## create Recur object for checking
    out <- methods::new("Recur",
                        .Data = sorted_mat[rev_ord, , drop = FALSE],
                        call = match.call(),
                        ID = ID,
                        ord = ord,
                        rev_ord = rev_ord,
                        first_idx = first_idx,
                        last_idx = last_idx,
                        check = check,
                        time_class = time_class)

    ## perform optional checks
    check_Recur(out, check = check)

    ## return
    out
}

##' Checks for Recurrent Event Data
##'
##' Perform several checks for recurrent event data and update object
##' attributions if some rows of the contained data (in the \code{.Data} slot)
##' have been removed by such as \code{na.action}.
##'
##' @inheritParams Recur
##' @param x An \code{Recur} object.
##'
##' @return An \code{Recur} object invisibly.
##' @export
check_Recur <- function(x, check = c("hard", "soft", "none"))
{
    ## is x of class Recur
    if (! is.Recur(x)) {
        stop("The 'x' must be an 'Recur' object.", call. = FALSE)
    }

    ## get check
    x@check <- check <- match.arg(check, choices = c("hard", "soft", "none"))

    ## determine whether some rows in x has been removed
    n_row <- nrow(x@.Data)
    is_obj_diff <- length(attr(x, "ord")) != n_row ||
        length(attr(x, "rev_ord")) != n_row

    ## if data input has an attr 'ID' for internal usage
    if (is_obj_diff) {
        ## check ID
        if (length(x@ID) > max(x[, "id"])) {
            stop("Found unknown ID.")
        }
        ## update ord and rev_ord
        x@ord <- order(x[, "id"], x[, "time2"], - x[, "event"])
        x@rev_ord <- order(x@ord)
        ## update first_idx and last_idx
        sorted_id <- x[x@ord, "id"]
        x@first_idx <- which(! duplicated.default(sorted_id))
        x@last_idx <- which(! duplicated.default(sorted_id, fromLast = TRUE))
        x@ID <- x@ID[sorted_id][x@first_idx]
    }
    ## early exit
    if (check == "none") {
        return(invisible(x))
    }

    ## sort data by id, time2, and - event
    sObj <- x@.Data[x@ord, , drop = FALSE]
    sID <- x@ID[sObj[, "id"]]
    sTime1 <- sObj[, "time1"]
    sTime2 <- sObj[, "time2"]
    sEvent <- sObj[, "event"]
    sTerminal <- sObj[, "terminal"]
    sCensor <- sEvent <= 0

    msg_fun <- if (check == "hard") { stop } else { warning }
    ## stop if event time after censoring time or without censoring time
    idx <- ! sCensor[x@last_idx]
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Subjects having events at or after censoring:",
            paste0(paste(sID[x@last_idx][idx], collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## stop if more than one terminal event time
    terminalID <- sID[sTerminal > 0]
    idx <- duplicated(terminalID)
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Subjects having multiple terminal events:",
            paste0(paste(terminalID[idx], collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## stop if any terminal event happens before the last time
    idx <- sTerminal[- x@last_idx] > 0
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Subjects having terminal events before censoring:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## the following check is not disabled for time-verying covariates
    ## stop if more than one censoring time
    ## cenID <- sID[sCensor]
    ## idx <- duplicated(cenID)
    ## if (any(idx)) {
    ##     msg_fun(wrapMessages(
    ##         "Subjects having multiple censoring times:",
    ##         paste0(paste(cenID[idx], collapse = ", "), ".")
    ##     ), call. = FALSE)
    ## }

    ## stop if missing value of 'time'
    idx <- is.na(sTime1) | is.na(sTime2)
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Missing times!",
            "Please check subject:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## 'time2' has to be later than the 'time1'
    idx <- sTime2 < sTime1
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Event times must be >= origin.",
            "Please check subject:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }

    ## 'time1' has to be not earlier than last 'time2'
    lag_sTime2 <- c(NA, sTime2[- n_row])
    lag_sTime2[x@first_idx] <- NA
    idx <-  ! is.na(lag_sTime2) & sTime1 < lag_sTime2
    if (any(idx)) {
        msg_fun(wrapMessages(
            "Recurrent episodes cannot be overlapped.",
            "Please check subject:",
            paste0(paste(unique(sID[idx]), collapse = ", "), ".")
        ), call. = FALSE)
    }
    ## return the (updated) x
    invisible(x)
}


##' Is the xect from the Recur class?
##'
##' Return \code{TRUE} if the specified xect is from the \code{\link{Recur}}
##' class, \code{FALSE} otherwise.
##'
##' @param x An \code{R} xect.
##' @return A logical value.
##' @export
is.Recur <- function(x)
{
    is(x, "Recur")
}


##' Convert An Recur Object to A Character Vector
##'
##' Summarize and convert the recurrent episodes for each subjects into
##' character strings.
##'
##' This function is intended to be a helper function for the `show()` method of
##' `Recur` objects.  To be precise, the function set the maximum number of
##' recurrent episodes for each subject to be `max(2L,
##' as.integer(getOption("reda.Recur.maxPrint")))`.  By default, at most three
##' recurrent episodes will be summarized for printing.  When subjects having
##' more than three recurrent episodes, the first
##' `getOption("reda.Recur.maxPrint") - 1` number of recurrent episodes and the
##' last one will be summarized.  One may use `options()` to adjust the setting.
##' For example, the default value is equivalent to `options(reda.Recur.maxPrint
##' = 3)`.
##'
##' @param x An Recur object.
##' @param ... Other arguments for future usage.
setMethod(f = "as.character", signature = "Recur",
          definition = function(x, ...) {
              ## determine the number of significant digits
              charNum <- unique(as.character(
                  x@.Data[, c("time1", "time2")]
              ))
              tmpList <- strsplit(charNum, "\\.")
              sigMax <- min(
                  max(sapply(tmpList, function(a) {
                      if (length(a) > 1)
                          return(nchar(a[2L]))
                      0
                  })),
                  max(3, getOption("digits") - 3)
              )
              fmt <- sprintf("(%s.%df, %s.%df%s]",
                             "%", sigMax, "%", sigMax, "%s")
              sorted_dat <- x@.Data[x@ord, , drop = FALSE]
              sorted_id <- x@ID[x@.Data[x@ord, "id"]]
              ## get options on max print
              max_print <- max(2L, as.integer(getOption("reda.Recur.maxPrint")))
              ## create a character vector representing the recurrent events
              char_rec <- tapply(
                  seq_along(x@ord), sorted_id,
                  function(idx) {
                      sub_time1 <- sorted_dat[idx, "time1"]
                      sub_time2 <- sorted_dat[idx, "time2"]
                      sub_is_censored <- sorted_dat[idx, "event"] == 0
                      sub_terminal <- max(sorted_dat[idx, "terminal"],
                                          na.rm = TRUE)
                      sub_end <- ifelse(sub_terminal > 0, "*", "+")
                      sub_sign <- ifelse(sub_is_censored, "+", "")
                      sub_sign[length(idx)] <- sub_end
                      out <- sprintf(fmt, sub_time1,
                                     sub_time2, sub_sign)
                      char_id <- sprintf("%s:", sorted_id[idx[1L]])
                      out_char <- if (length(sub_time1) > max_print) {
                                      paste(c(out[seq_len(max_print - 1)],
                                              "...",
                                              out[length(out)]),
                                              collapse = ", ")
                                  } else {
                                      paste(out, collapse = ", ")
                                  }
                      paste(char_id, out_char)
                  })
              unname(as.character(char_rec[unique(sorted_id)]))
          })


## helper function to process 'time'
process_time <- function(x, allow_null = TRUE) {
    ## skip if x is null or a list
    if (is.null(x)) {
        if (allow_null) {
            return(x)
        } # else
        stop("The 'time' cannot be 'NULL'.", call. = FALSE)
    }
    if (is.list(x))
        return(x)
    if (inherits(x, "difftime") || inherits(x, "Date") || inherits(x, "POSIXt"))
        x <- as.numeric(x)
    if (! is.numeric(x))
        stop("Invalid 'time' variable.  See '?Recur' for details.",
             call. = FALSE)
    ## return
    x
}
