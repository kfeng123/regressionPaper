library(MASS)
library(CVXR)
library(SILM)
library(scalreg)


set.seed(1)


lassoTest_on <- TRUE
lassoTestOracle_on <- TRUE

RepTime <- 200
M <- 2000
alpha <- 0.05
plotPointSize <- 100

q <- 10
p <- 100

SNRsequence <- c(0,10,20,30)

betaDistribution <- "unif"


XGen_list <- c("uniform","equalCor","factor", "rowCor")
n_list <- c(50)
epsilonDis_list <- c("t","chi")
betabGen_list <- c("dense","sparse")


for(XGen in XGen_list)
    for(n in n_list)
        for(epsilonDis in epsilonDis_list)
            for(betabGen in betabGen_list){
                source("./code.R")
                source("./plot.R")
            }



# for(XGen in c("rowCor"))
#     for(n in c(50))
#         for(epsilonDis in c("t"))
#             for(betabGen in c("dense"))
#                 source("./code.R")
# 
# 
# for(XGen in c("rowCor"))
#     for(n in c(50))
#         for(epsilonDis in c("t"))
#             for(betabGen in c("dense")){
#                 source("./plot.R")
#             }

