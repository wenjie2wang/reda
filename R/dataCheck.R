### internal function ==========================================================
## collation after class.R
#' @include class.R
check_Survr <- function(dat) {
    ## check missing value on 'ID'
    if (any(is.na(dat$ID))) {
        stop("'ID' cannot be missing.")
    }
    ## check coding and missing value on 'event'
    if (any(! dat$event %in% 0:1)) {
        stop("'event' must be coded as 0 (censoring) or 1 (event).")
    }
    ## if dat input has an attr 'ID'
    nID <- attr(dat, "ID")
    if (! is.null(nID)) {
        dat$IDnam <- nID
    } else {
        ## check whether 'ID' is numeric or not. convert if not.
        dat$IDnam <- factor(dat$ID, levels = unique(dat$ID))
        dat$ID <- as.numeric(dat$IDnam)
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
        } else {2}
        mis_time0 <- if (length(time0) > 0) {
            ## missing indicator of time for event == 0
            ifelse(any(is.na(time0)), 1, 0)
        } else {2}
        ## issue #1: without censoring time or more than one censoring time
        censor1 <- ifelse(sum(subdat$event == 0, na.rm = TRUE) != 1, 1, 0)
        ## issue #2: event time after censoring time
        censor2 <- if (mis_time1 == 0 && mis_time0 == 0) {
            ifelse(max(time1) >= min(time0), 1, 0)
        } else {2}
        ## return
        cbind(subdat, mis_time1, mis_time0, censor1, censor2)
    }
    outDat <- plyr::ddply(dat, "ID", check_ddply)
    ## stop if missing value of 'time' for event == 1
    ID_mis_time1 <- with(subset(outDat, mis_time1 == 1), unique(IDnam))
    if (length(ID_mis_time1) > 0) {
        stop(paste("There is missing value on event time for subject:", 
                   paste0(ID_mis_time1, collapse = ", ")))
    }
    ## stop if missing value of 'time' for event == 0
    ID_mis_time0 <- with(subset(outDat, mis_time0 == 1), unique(IDnam))
    if (length(ID_mis_time0) > 0) {
        stop(paste("Censoring time is missing for subject:", 
                   paste0(ID_mis_time0, collapse = ", ")))
    }
    ## stop if no censoring time or more than one censoring time
    ID_censor1 <- with(subset(outDat, censor1 == 1), unique(IDnam))
    if (length(ID_censor1) > 0) {
        message("Every subject must have one (and only one) censored time.")
        stop(paste("Check subject: ",
            paste0(ID_censor1, collapse = ", ")))
    }
    ## stop if event time after censoring time
    ID_censor2 <- with(subset(outDat, censor2 == 1), unique(IDnam))
    if (length(ID_censor2) > 0) {
        message("Event time should be earlier than censoring time.")
        stop(paste("Check subject:",
                   paste0(ID_censor2, collapse = ", ")))
    }
    ## return
    out <- outDat[, c("ID", "time", "event")]
    attr(out, "ID") <- outDat$IDnam
    invisible(out)
}
