library(reda)

data(simuDat)
simuDat$CID <- paste0("A", simuDat$ID)
simuDat$origin <- 0

## error if ID, time, or event is missing
expect_error(with(simuDat, Survr()), "ID")
expect_error(with(simuDat, Survr(CID, )), "Time")
expect_error(with(simuDat, Survr(ID, time)), "Event")

## error if ID contains missing values
tmpDat <- simuDat
tmpDat[500, "CID"] <- NA
expect_error(with(tmpDat, Survr(CID, time, event)),
             "missing values")

## error if time or origin cannot be converted to numerical values
expect_error(with(simuDat, Survr(CID, as.character(time), event)),
             "Time")
expect_error(with(simuDat, Survr(CID, time, event, origin = "0")),
             "Origin")

## issue 1: event time after censoring time or without censoring time
tmpDat <- simuDat
tmpDat[c(5, 10), "event"] <- 1
expect_error(with(tmpDat, Survr(CID, time, event)),
             "A1, A4.")

tmpDat <- simuDat
tmpDat[c(27, 30), "time"] <- 1
expect_error(with(tmpDat, Survr(CID, time, event)),
             "A5, A6.")

## issue 2: more than one censoring time
tmpDat <- simuDat
tmpDat[c(25, 28), "event"] <- 0
expect_error(with(tmpDat, Survr(CID, time, event)),
             "A5, A6.")

## error if time contains missing values
tmpDat <- simuDat
tmpDat[500, "time"] <- NA_real_
expect_error(with(tmpDat, Survr(CID, time, event)),
             "A100.")

## error if origin contains missing values
tmpDat <- simuDat
tmpDat[490, "origin"] <- NA_real_
expect_error(with(tmpDat, Survr(CID, time, event, origin)),
             "A99.")

## error if time is earlier than origin time
tmpDat <- simuDat
tmpDat[490, "origin"] <- 200
expect_error(with(tmpDat, Survr(CID, time, event, origin)),
             "A99.")

## error if one subject has different origins
tmpDat <- simuDat
tmpDat[489 : 491, "origin"] <- - 1
expect_error(with(tmpDat, Survr(CID, time, event, origin)),
             "A100, A98.")

## error caused after moving observation containing missing values
tmpDat <- simuDat
tmpDat[c(8, 94), "x1"] <- NA
expect_error(rateReg(Survr(CID, time, event) ~ x1, tmpDat,
                     control = list(verbose = FALSE)),
             "A16, A3.")
