library(reda)

set.seed(123)
data(simuDat)
simuDat$CID <- paste0("A", simuDat$ID)
simuDat <- simuDat[sample(nrow(simuDat), nrow(simuDat)), ]

## when only `time` is given
ex1 <- Recur(3:5, check = "soft")
expect_equal(ex1[, "time1"], rep(0, 3))
expect_equal(ex1[, "time2"], 3:5)
expect_equal(ex1[, "id"], 1:3)
expect_equal(ex1[, "event"], rep(0, 3))
expect_equal(ex1[, "terminal"], rep(0, 3))
expect_equal(ex1[, "origin"], rep(0, 3))
expect_equal(ex1@first_idx, 1:3)
expect_equal(ex1@last_idx, 1:3)
expect_equal(ex1@ord, 1:3)
expect_equal(ex1@rev_ord, 1:3)
expect_equal(ex1@check, "soft")

## When `time` and `id` are given
ex2 <- Recur(6:1, id = rep(1:2, 3), check = "none")
expect_equal(ex2[, "time1"], c(4, 3, 2, 1, 0, 0))
expect_equal(ex2[, "time2"], 6:1)
expect_equal(ex2[, "id"], rep(1:2, 3))
expect_equal(ex2[, "event"], c(rep(0, 2), rep(1, 4)))
expect_equal(ex2[, "terminal"], rep(0, 6))
expect_equal(ex2[, "origin"], rep(0, 6))
expect_equal(ex2@first_idx, c(1, 4))
expect_equal(ex2@last_idx, c(3, 6))
expect_equal(ex2@ord, c(5, 3, 1, 6, 4, 2))
expect_equal(ex2@rev_ord, c(3, 6, 2, 5, 1, 4))
expect_equal(ex2@check, "none")

## Helper `%to%` (or `%2%`) for `time`
time1 <- c(1, 5, 7)
time2 <- c(3, 7, 9)
my_id <- c("A1", "A1", "A2")
ex3 <- Recur(time1 %to% time2, id = my_id)
expect_equal(ex3[, "time1"], time1)
expect_equal(ex3[, "time2"], time2)
expect_equal(ex3[, "id"], as.numeric(factor(my_id)))
expect_equal(ex3[, "event"], c(1, 0, 0))
expect_equal(ex3[, "terminal"], rep(0, length(time1)))
expect_equal(ex3[, "origin"], c(1, 1, 7))
expect_equal(ex3@first_idx, c(1, 3))
expect_equal(ex3@last_idx, c(2, 3))
expect_equal(ex3@ord, 1:3)
expect_equal(ex3@rev_ord, 1:3)
expect_equal(ex3@check, "hard")

ex4 <- Recur(list(time1 = time1, time2 = time2), id = c("A1", "A1", "A2"))
expect_equivalent(ex3, ex4)

## with `origin` and `terminal`
ex5 <- Recur(3:5, origin = 1, terminal = 1)
expect_equal(ex5[, "time1"], rep(1, nrow(ex5)))
expect_equal(ex5[, "time2"], 3:5)
expect_equal(ex5[, "id"], 1:3)
expect_equal(ex5[, "event"], rep(0, nrow(ex5)))
expect_equal(ex5[, "terminal"], rep(1, nrow(ex5)))
expect_equal(ex5[, "origin"], rep(1, nrow(ex5)))
expect_equal(ex5@first_idx, seq_len(nrow(ex5)))
expect_equal(ex5@last_idx, seq_len(nrow(ex5)))
expect_equal(ex5@ord, seq_len(nrow(ex5)))
expect_equal(ex5@rev_ord, seq_len(nrow(ex5)))
expect_equal(ex5@check, "hard")

ex6 <- Recur(3:5, id = c("A1", "A1", "A2"), origin = 1:2, terminal = c(0, 1))
ex7 <- Recur(3:5, id = c("A1", "A1", "A2"),
             origin = c(1, 1, 2), terminal = c(0, 0, 1))
expect_equivalent(ex6, ex7)

## error due to inappropriate length
expect_error(Recur(1:10, origin = c(1, 2)))
expect_error(Recur(1:10, terminal = c(1, 2)))

## data checking rules

## Every subject must have one censoring not before any event time
expect_error(Recur(1:3, id = rep("A1", 3), event = c(0, 0, 1)), "A1")

## Every subject must have one terminal event time.
expect_error(Recur(1:3, id = rep("A1", 3), terminal = c(0, 1, 1)), "A1")

## Recurrent episode without events is allowed for time-varying covariates
ex8 <- Recur(1:3, id = rep("A1", 3), event = c(0, 1, 0))
expect_equal(ex8[, "time1"], c(0, 1, 2))
expect_equal(ex8[, "time2"], 1:3)
expect_equal(ex8[, "id"], rep(1, nrow(ex8)))
expect_equal(ex8[, "event"], c(0, 1, 0))
expect_equal(ex8[, "terminal"], rep(0, nrow(ex8)))
expect_equal(ex8@ord, 1:3)
expect_equal(ex8@rev_ord, 1:3)
expect_equal(ex8@first_idx, 1)
expect_equal(ex8@last_idx, 3)
expect_equal(ex8@check, "hard")

## Event or censoring times cannot be missing
expect_error(Recur(c(1:2, NA), id = rep("A1", 3)))

## Event times cannot be earlier than the origin time
expect_error(Recur(3:5, id = rep("A1", 3), origin = 10))
expect_error(Recur(3:5 %to% 1:3, id = rep("A1", 3)))

## The recurrent episode cannot be overlapped
expect_error(Recur(c(0, 3, 5) %to% c(1, 6, 10), id = rep("A1", 3)))
