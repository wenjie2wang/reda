testdat <- data.frame(matrix(c(1, 20, 1, 5, 1, 70, 1, 5, 1, 100, 0, 
                    5, 2, 60, 1, 8, 2, 100, 0, 8), nrow = 5, byrow = TRUE))
colnames(testdat) <- c("ID", "Time", "Event", "X")
baselinepieces <- c(50, 100)
ini <- c(0.3, 0.5, 0.15, 0.16)
testdat
#   ID Time Event X
# 1  1   20     1 5
# 2  1   70     1 5
# 3  1  100     0 5
# 4  2   60     1 8
# 5  2  100     0 8


## dataset: colon from package survrec
require(survrec)
data(MMC)
testdat <- MMC
colnames(testdat)[1:3] <- c("ID", "Time", "Event")
testdat$group <- as.numeric(testdat$group)
BaselinePieces <- quantile(testdat$Time, probs = c(0.25, 0.5, 0.75))
attr(BaselinePieces, "names") <- NULL
BaselinePieces <- c(BaselinePieces, max(testdat$Time))
BaselinePieces
# [1]  53  85 121 284
ini <- rep(0.15, 4)
