testdat <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 2), rep(4, 3)),
                      Time = c(20, 70, 100, 60, 100, 50, 70, 30, 50, 75),
                      Event = c(1, 1, 0, 1, 0, 1, 0, 1, 1, 0), 
                      X1 = c(rep(5, 3), rep(8, 2), rep(7, 2), rep(NA, 3)),
                      X2 = gl(n = 2, k = 5, length = 10))
testdat
#    ID Time Event X1 X2
# 1   1   20     1  5  1
# 2   1   70     1  5  1
# 3   1  100     0  5  1
# 4   2   60     1  8  1
# 5   2  100     0  8  1
# 6   3   50     1  7  2
# 7   3   70     0  7  2
# 8   4   30     1 NA  2
# 9   4   50     1 NA  2
# 10  4   75     0 NA  2

## test
heartfit <- heart(formula = survrec::Survr(ID, Time, Event) ~ X1 + X2, 
                  data = testdat, subset = ID %in% 1:4, 
                  na.action = na.omit, 
                  contrasts = list("contr.sum"),
                  baselinepieces = c(50, 100), 
                  start = list(beta = c(0.3, 1), 
                               theta = 0.5, 
                               alpha = c(0.15, 0.16)),
                  control = list(gradtol = 1e-6, stepmax = 1e5, 
                                 steptol = 1e-6, iterlim = 1e2))
show(heartfit) # or simply call: heartfit
str(heartfit)


## dataset: colon from package survrec
require(survrec)
data(colon)
testdat <- colon
BaselinePieces <- quantile(testdat$time, probs = c(0.25, 0.5, 0.75, 1))
attr(BaselinePieces, "names") <- NULL
BaselinePieces
# [1]  21  216  880 2175




