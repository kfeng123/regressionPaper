set.seed(1)


RepTime <- 2500
M <- 800
alpha <- 0.05
plotPointSize <- 100

n <- 100
q <- 10
p <- 1000


betaDistribution <- "unif"
XGen <- "equalCor"
epsilonDis <- "t"
betabGen <- "dense"

for(XGen in c("factor","equalCor","uniform"))
    for(n in c(50,100))
        for(epsilonDis in c("t","chi"))
            for(betabGen in c("dense","sparse"))
                source("./code.R")



