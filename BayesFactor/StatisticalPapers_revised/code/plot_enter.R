library(ggplot2)

for(XGen in c("uniform","equalCor","factor"))
    for(n in c(50,100))
        for(epsilonDis in c("t","chi"))
            for(betabGen in c("dense","sparse")){
                source("./plot.R")
            }