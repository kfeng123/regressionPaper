library(MASS)
library(CVXR)
library(SILM)
library(scalreg)


set.seed(1)


lassoTest_on <- FALSE

RepTime <- 3000
M <- 2000
alpha <- 0.05
plotPointSize <- 100

q <- 10
p <- 1000

SNRsequence <- c(0,10,20,30)

betaDistribution <- "unif"
XGen <- "equalCor"
epsilonDis <- "t"
betabGen <- "dense"

for(XGen in c("uniform","equalCor","factor"))
    for(n in c(50, 100))
        for(epsilonDis in c("t","chi"))
            for(betabGen in c("dense"))
                source("./code.R")


