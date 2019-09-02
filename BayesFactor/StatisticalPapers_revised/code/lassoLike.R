n <- 100
p <- 1000

X <- rnorm(n*p); dim(X) <- c(n,p)
beta <- rep(1,p)
#beta[sample(p,p-5, replace = TRUE)] <- 0
y <- rnorm(n,0,2)  
############ Code provided by Xianyang Zhang

n <- 50
p <- 300

X <- rnorm(n*p); dim(X) <- c(n,p)
beta <- rep(1,p)
#beta[sample(p,p-5, replace = TRUE)] <- 0
set <- 1:5
source("./ZhangXianYang/covarianceMatrix.R")
lassoTest <- rep(0,1000)
for(i in 1:1000){
    y <- rnorm(n,0,2)  
    source("./ZhangXianYang/lassoTest.R")
    lassoTest[i] <- lassoTest_one_result
}
