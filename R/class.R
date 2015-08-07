################################################################################
##
##   R package heart by Haoda Fu, Jun Yan, and Wenjie Wang
##   Copyright (C) 2015
##
##   This file is part of the R package heart.
##
##   The R package heart is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package heart is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package heart. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################


## create S3 Class called "survr" as formula response for recurrent event data
#' An S3 Class to Represent Formula Response for Recurrent Event Data
#' 
#' \code{Survr} is an S3 class to represent 
#' formula response for recurrent event data. 
#' @param ID numeric identificator of each subject. 
#' For character 'ID', please use \code{\link[stats]{setNames}} to 
#' create a named numeric vector. For example, if original 'ID' is 
#' \code{c("A", "B", "C")}, use 
#' \code{ID <- setNames(1:3, LETTERS[1:3])} as new numeric 'ID'.
#' @param time time of reccurence. For each subject the last time are censored.
#' @param event the status indicator, 
#' 0 = censored, 1 = recurrent event. 
#' @aliases Survr
#' @seealso \code{\link{heart}}
#' @importFrom plyr ddply
#' @export
Survr <- function (ID, time, event) {
  dat <- data.frame(ID, time, event)
  check_Survr(dat)
  ans <- as.matrix(cbind(ID, time, event))
  oldClass(ans) <- "Survr"
  invisible(ans)
}


## create S4 Class called "heart" for heart object from function heart
#' An S4 Class to Represent a Fitted HEART Model
#' 
#' \code{heart-class} is an S4 class to represent a HEART model fits. 
#' \code{\link{heart}} produces objects of this class.  
#' See ``Slots'' for details.
#' 
#' @slot call function call.
#' @slot formula formula.
#' @slot baselinepieces a numeric vector.
#' @slot estimates list.
#' @slot control list.
#' @slot start list.
#' @slot na.action a length-one character vector.
#' @slot xlevels list.
#' @slot contrasts list.
#' @slot convergence an integer.
#' @slot hessian a numeric matrix.
#' @aliases heart-class
#' @seealso \code{\link{heart}}
#' @export
setClass(Class = "heart", 
         slots = c(call = "call", 
                   formula = "formula", 
                   baselinepieces = "numeric",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergence = "integer", 
                   hessian = "matrix"))


## create S4 Class called "summary.heart" for summary.heart object from summary
#' An S4 Class to Represent Summary of heart-class Object
#' 
#' summary.heart-class is an S4 class with selective slots 
#' of heart-class object.  See ``Slots'' for details.  
#' \code{\link{summary}} produces objects of this class. 
#'  
#' @slot call function call.
#' @slot baselinepieces a numeric vector.
#' @slot coefficients a numeric matrix.
#' @slot theta numeric a matrix.
#' @slot baseline a numeric matrix.
#' @aliases summary.heart-class
#' @seealso \code{\link{summary,heart-method}} 
#' @export
setClass(Class = "summary.heart", 
         slots = c(call = "call", 
                   baselinepieces = "numeric",
                   coefficients = "matrix",
                   theta = "matrix",
                   baseline = "matrix"))


#' An S4 Class to Represent Computed Empirical MCF
#' 
#' An S4 class to represent computed empirical mean cumulative function (MCF).
#' \code{\link{MCF}} produces objects of this class.  
#' @slot call function call
#' @slot formula formula. 
#' @slot MCF a data.frame.
#' @slot multigroup a logical value.
#' @aliases empirMCF-class
#' @seealso \code{\link{MCF,formula-method}}
#' @importFrom methods setClass
#' @export
setClass(Class = "empirMCF", 
         slots = c(call = "call", formula = "formula", MCF = "data.frame", 
                   multigroup = "logical"))


#' An S4 Class to Represent Estimated MCF from HEART Model
#' 
#' An S4 class to represent estimated mean cumulative function (MCF) 
#' from HEART Model.
#' \code{\link{MCF}} produces objects of this class.  
#' 
#' @slot formula formula.
#' @slot baselinepieces a numeric vector.
#' @slot newdata a numeric matrix.
#' @slot MCF a data.frame.
#' @slot level a numeric value between 0 and 1.
#' @slot na.action a length-one character vector.
#' @slot control list.
#' @slot multigroup a logical value. 
#' @aliases heartMCF-class
#' @seealso \code{\link{MCF,heart-method}}
#' @export
setClass(Class = "heartMCF", 
         slots = c(formula = "formula", baselinepieces = "numeric",
                   newdata = "matrix", MCF = "data.frame", level = "numeric", 
                   na.action = "character", control = "list", 
                   multigroup = "logical"))


### internal function ========================================================== 
check_Survr <- function(dat) {
  ## check missing value on 'ID'
  if (any(is.na(dat$ID))) {
    stop("'ID' cannot be missing.")
  }
  ## check coding and missing value on 'event'
  if (any(! dat$event %in% 0:1)) {
    stop("'event' must be coded as 0 (censoring) or 1 (event).")
  }
  ## nonsense, just to suppress Note from R CMD check --as-cran
  mis_time1 <- mis_time0 <- censor1 <- censor2 <- event <- NULL
  ## check function
  check_ddply <- function (subdat) {
    subdat <- subdat[order(subdat$time), ]
    ## check missing values on 'time'
    time1 <- with(subset(subdat, event == 1), time)
    time0 <- with(subset(subdat, event == 0), time)
    mis_time1 <- if (length(time1) > 0) {
      ## missing indicator of time for event == 1
      ifelse(any(is.na(time1)), 1, 0)
    } else { 2 }
    mis_time0 <- if (length(time0) > 0) {
      ## missing indicator of time for event == 0
      ifelse(any(is.na(time0)), 1, 0)
    } else { 2 }
    ## issue #1: without censoring time or more than one censoring time
    censor1 <- ifelse(sum(subdat$event == 0, na.rm = TRUE) != 1, 1, 0)
    ## issue #2: event time after censoring time
    censor2 <- if (mis_time1 == 0 && mis_time0 == 0) {
      ifelse(max(time1) >= min(time0), 1, 0)
    } else { 2 }
    ## return
    cbind(subdat, mis_time1, mis_time0, censor1, censor2)
  }
  outdat <- plyr::ddply(dat, "ID", check_ddply)
  ## stop if missing value of 'time' for event == 1
  ID_mis_time1 <- with(subset(outdat, mis_time1 == 1), unique(ID))
  if (length(ID_mis_time1) > 0) {
    stop(paste("There is missing value on event time for subject:", 
               paste0(ID_mis_time1, collapse = ", ")))
  }
  ## stop if missing value of 'time' for event == 0
  ID_mis_time0 <- with(subset(outdat, mis_time0 == 1), unique(ID))
  if (length(ID_mis_time0) > 0) {
    stop(paste("Censoring time is missing for subject:", 
               paste0(ID_mis_time0, collapse = ", ")))
  }
  ## stop if no censoring time or more than one censoring time
  ID_censor1 <- with(subset(outdat, censor1 == 1), unique(ID))
  if (length(ID_censor1) > 0) {
    stop(paste(
      "Every subject must have one (and only one) censored time. Check subject:",
      paste0(ID_censor1, collapse = ", ")))
  }
  ## stop if event time after censoring time
  ID_censor2 <- with(subset(outdat, censor2 == 1), unique(ID))
  if (length(ID_censor2) > 0) {
    stop(paste(
      "Event time should be earlier than censoring time. Check subject:",
      paste0(ID_censor2, collapse = ", ")))
  }
  ## return
  return(NULL)
}
