testdat <- matrix(c(1, 20, 1, 5, 1, 70, 1, 5, 1, 100, 0, 
                    5, 2, 60, 1, 8, 2, 100, 0, 8), nrow = 5, byrow = TRUE)
colnames(testdat) <- c("ID", "Time", "Event", "X")
BLpieces <- c(50, 100)
testdat
#      ID Time Event X
# [1,]  1   20     1 5
# [2,]  1   70     1 5
# [3,]  1  100     0 5
# [4,]  2   60     1 8
# [5,]  2  100     0 8
