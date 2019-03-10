

RepTime <- 50
M <- 2000
alpha <- 0.05
plotPointSize <- 100

n <- 100
q <- 10
p <- 1000


betaDistribution <- "unif"
XGen <- "equalCor"
epsilonDis <- "t"
betabGen <- "dense"

for(XGen in c("uniform","equalCor","factor"))
    for(n in c(100))
        for(epsilonDis in c("t","chi"))
            for(betabGen in c("dense","sparse"))
                source("./code.R")



