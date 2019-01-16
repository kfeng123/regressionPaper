library(xtable)
set.seed(1)

RepTime <- 100
M <- 1000
alpha <- 0.05

n <- 100
q <- 10
p <- 1000


XGen <- "equalCor"

n <- 100

epsilonDis <- "t"

betabGen <- "dense"
source("./code.R")
biaoge1 <-as.matrix(ohResult)

betabGen <- "sparse"
source("./code.R")
biaoge2 <- as.matrix(ohResult)

epsilonDis <- "chi"

betabGen <- "dense"
source("./code.R")
biaoge3 <-as.matrix(ohResult)

betabGen <- "sparse"
source("./code.R")
biaoge4 <- as.matrix(ohResult)


finalBiaoge <- cbind(biaoge1,biaoge2[,-1],biaoge3[,-1],biaoge4[,-1])
print(xtable(finalBiaoge,auto=TRUE),include.rownames = FALSE)



