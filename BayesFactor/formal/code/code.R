library(MASS)

n <- 50
q <- 10
p <- 100
Xa <- rnorm(n*q)
dim(Xa) <- c(n,q)
Xb <- rnorm(n*p)
dim(Xb) <- c(n,p)
y <- rep(1,n)


proposedTest <- function(y,Xa,Xb){
    n <- nrow(Xa)
    q <- ncol(Xa)
    p <- ncol(Xb)
    tildeUa <- MASS::Null(Xa)
    XbXbT <- Xb %*% t(Xb)
    XbStarXbStarT <- t(tildeUa) %*% XbXbT %*% tildeUa
    tmp1 <- tildeUa %*% solve(XbStarXbStarT) %*% t(tildeUa)
    theNumerator <- as.numeric( 
        t(y) %*%  %*% y
        )
    theDenominator <- sum(
                            (t(tildeUa) %*% y)^2
                          )
    proposedStat <- - theNumerator / theDenominator
     
    
}