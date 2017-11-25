context("Testing mcf,formula method")

test_that("test mcf,formula method", {
    data(valveSeats)
    valveSeats$group_ <- cut(valveSeats$ID, c(250, 400, 450))

    ## test na.action
    testDat <- valveSeats
    testDat[1, "group_"] <- NA
    test_mcf <- mcf(Survr(ID, Days, No.) ~ group_,
                    testDat, na.action = "na.exclude",
                    variance = "bootstrap",
                    control = list(B = 30))
    expect_equal(nrow(valveSeats) - nrow(test_mcf@data), 1L)

    ## test bootstrap variance
    expect_equivalent(class(
        mcf(Survr(ID, Days, No.) ~ group_, valveSeats,
            variance = "bootstrap", logConfInt = TRUE,
            control = list(B = 30, se.method = "normality"))
    ), "mcf.formula")
    expect_equivalent(class(
        mcf(Survr(ID, Days, No.) ~ group_, valveSeats,
            variance = "bootstrap",
            control = list(B = 30, se.method = "normality",
                           ci.method = "percentile"))
    ), "mcf.formula")
    expect_equivalent(class(
        mcf(Survr(ID, Days, No.) ~ group_, valveSeats,
            variance = "bootstrap", logConfInt = TRUE,
            control = list(B = 30, se.method = "normality",
                           ci.method = "percentile"))
    ), "mcf.formula")

    ## test its plot-method
    mcf0 <- mcf(Survr(ID, Days, No.) ~ 1, valveSeats)
    expect_equivalent(class(
        plot(mcf0, conf.int = TRUE, addOrigin = TRUE, mark.time = TRUE)
    ), c("gg", "ggplot"))
    expect_equivalent(class(
        plot(test_mcf, conf.int = TRUE, addOrigin = TRUE, mark.time = TRUE,
             legendLevels = c("G1", "G2"), legendName = "Group")
    ), c("gg", "ggplot"))
    expect_equivalent(class(
        plot(test_mcf, lty = 1:2, col = 1:2)
    ), c("gg", "ggplot"))

    ## show method
    expect_output(show(mcf0), "MCF")
})
