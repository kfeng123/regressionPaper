library(glmnet)
library(scalreg)
n <- 100
p <- 1000

X <- rnorm(n*p); dim(X) <- c(n,p)
beta <- rep(1,p)
#beta[sample(p,p-5, replace = TRUE)] <- 0
y <- rnorm(n,0,2)  
haha <- scalreg(X,y)

lasso_model <- glmnet(X, y,
                        family="gaussian",
                        alpha=1,
                        nlambda = 1000,
                        standardize= FALSE,
                        intercept=FALSE)
theCoef <- coef(lasso_model, s = 0.1)
thePred <- predict(lasso_model, X, s=0.1)


############ Scaled LASSO ########################
SL.iterTime <- 1000
SL.lambda_array <- rep(0, SL.iterTime)
SL.lambda_array[1] <- 0.00001
SL.sigma_array <- rep(0, SL.iterTime)
for(i in 2:SL.iterTime){
    SL.sigma_est <- sqrt( mean( (y - predict(lasso_model, X, s= SL.lambda_array[i-1]) )^2 ) )
    SL.lambda_array[i] <- SL.lambda_array[i-1] * SL.sigma_est
    SL.sigma_array[i] <- SL.sigma_est
}



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
