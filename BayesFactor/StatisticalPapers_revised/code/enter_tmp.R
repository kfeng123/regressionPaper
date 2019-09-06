library(MASS)
library(CVXR)
library(SILM)
library(scalreg)


set.seed(1)


lassoTest_on <- TRUE

RepTime <- 200
M <- 2000
alpha <- 0.05
plotPointSize <- 100

n <- 50
q <- 10
p <- 1000

SNRsequence <- c(0,20)

betaDistribution <- "unif"

for(XGen in c("rowCor"))
    for(n in c(50))
        for(epsilonDis in c("t"))
            for(betabGen in c("dense"))
                source("./code.R")


for(XGen in c("rowCor"))
    for(n in c(50))
        for(epsilonDis in c("t"))
            for(betabGen in c("dense")){
                source("./plot.R")
            }

