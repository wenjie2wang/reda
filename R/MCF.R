## sample mean cumulative function
## return a data frame
#' @export
sample_MCF <- function(formula, data, subset, na.action) {
  ## record the function call 
  Call <- match.call()
  ## arguments check
  if(missing(formula)) {
    stop("formula argument is missing.")
  } 
  if(missing(data)) {
    data <- environment(formula)
  }
  if (! with(data, inherits(eval(Call[[2]][[2]]), "Survr"))) {
    stop("formula must be a survival recurrent object.")
  }
  ## Prepare data: ID, time, event ~ X
  mcall <- match.call(expand.dots = FALSE)
  mmcall <- match(c("formula", "data", "subset", "na.action"), names(mcall), 0L)
  mcall <- mcall[c(1L, mmcall)]
  ## drop unused levels in factors 
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- quote(stats::model.frame)
  mf <- eval(mcall, parent.frame())
  mm <- stats::model.matrix(formula, data = mf)
  ## number of covariates excluding intercept
  nbeta <- ncol(mm) - 1 
  ## covariates' names
  covar_names <- base::colnames(mm)[-1]
  ## data 
  dat <- as.data.frame(cbind(mf[, 1][, 1:3], mm[, -1]))
  colnames(dat) <- c("ID", "time", "event", covar_names)
  ## compute sample MCF 
  if (nbeta == 0) {
    num_pat <- length(unique(dat$ID))
    unitimes <- sort(unique(dat$time), decreasing = FALSE)
    unitimes <- head(unitimes, length(unitimes) - 1)
    num_at_risk <- apply(array(unitimes), 1, 
                         function(a) {
                           with(base::subset(dat, time <= a), 
                                num_pat - sum(event == 0))
                         })
    num_fails <- apply(array(unitimes), 1, 
                       function(a) {
                         with(base::subset(dat, time <= a), sum(event == 1))
                       })
    increment <- num_fails / num_at_risk
    sMCF <- cumsum(increment)
    ## data.frame to return
    num_at_risk <- c(length(unique(dat$ID)), num_at_risk)
    num_fails <- c(0, num_fails)
    increment <- c(0, increment)
    sMCF <- c(0, sMCF)
    unitimes <- c(0, unitimes)
    return(data.frame(time = unitimes, fails = num_fails, risk = num_at_risk, 
                      increment = increment, sMCF = sMCF))
  } else {
    ## argument check
    if (nbeta != 1) {
      stop("The covariate in formula must be a factor or '1'")
    } 
    ## extract the covariate
    fac <- with(data, eval(Call[[2]][[3]]))
    if (! is.factor(fac)){
      stop("The covariate in formula must be a factor or '1'")    
    }
    ## number of levels
    num_levels <- length(levels(fac))
    outdat <- NULL
    for (i in seq(num_levels)) {
      subdat <- base::subset(dat, subset = fac %in% levels(fac)[i])
      num_pat <- length(unique(subdat$ID))
      unitimes <- sort(unique(subdat$time), decreasing = FALSE)
      unitimes <- head(unitimes, length(unitimes) - 1)
      num_at_risk <- apply(array(unitimes), 1, 
                           function(a) {
                             with(base::subset(subdat, time <= a), 
                                  num_pat - sum(event == 0))
                           })
      num_fails <- apply(array(unitimes), 1, 
                         function(a) {
                           with(base::subset(subdat, time <= a), sum(event == 1))
                         })
      increment <- num_fails / num_at_risk
      sMCF <- cumsum(increment)
      ## subdata.frame to return
      num_at_risk <- c(length(unique(subdat$ID)), num_at_risk)
      num_fails <- c(0, num_fails)
      increment <- c(0, increment)
      sMCF <- c(0, sMCF)
      unitimes <- c(0, unitimes)
      outdat <- rbind(outdat, data.frame(unitimes, num_fails, num_at_risk, 
                                         increment, sMCF, levels(fac)[i]))
    }
    colnames(outdat) <- c("time", "fails", "risk", 
                          "increment", "sMCF", as.character(Call[[2]][[3]]))
    return(outdat)
  }
}
