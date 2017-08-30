context("Testing data checking procedure for Survr")

test_that("call reda::Survr", {
    data(simuDat)
    simuDat$CID <- paste0("A", simuDat$ID)

    ## issue 1: event time after censoring time or without censoring time
    tmpDat <- simuDat
    tmpDat[c(5, 10), "event"] <- 1
    expect_error(with(tmpDat, Survr(CID, time, event)),
                 "Please check subject: A1, A4.", fixed = TRUE)

    tmpDat <- simuDat
    tmpDat[c(27, 30), "time"] <- 1
    expect_error(with(tmpDat, Survr(CID, time, event)),
                 "Please check subject: A5, A6.", fixed = TRUE)

    ## issue 2: more than one censoring time
    tmpDat <- simuDat
    tmpDat[c(25, 28), "event"] <- 0
    expect_error(with(tmpDat, Survr(CID, time, event)),
                 "Please check subject: A5, A6.", fixed = TRUE)

    ## error caused after moving observation containing missing values
    tmpDat <- simuDat
    tmpDat[c(8, 94), "x1"] <- NA
    expect_error(rateReg(Survr(CID, time, event) ~ x1, tmpDat),
                 "Please check subject: A16, A3.", fixed = TRUE)

})
