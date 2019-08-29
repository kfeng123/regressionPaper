library(glmnet)
n <- 100
p <- 1000

X <- rnorm(n*p); dim(X) <- c(n,p)
y <- rnorm(n)
lasso_compute <- glmnet(X, y,
                        family="gaussian",
                        alpha=1,
                        nlambda = 1000,
                        standardize= FALSE,
                        intercept=FALSE)
theCoef <- coef(lasso_compute, s = 0.1)
thePred <- predict(lasso_compute, X, s=0.1)
