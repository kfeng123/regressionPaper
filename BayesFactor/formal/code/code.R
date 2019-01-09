library(MASS)

n <- 50
q <- 10
p <- 200
Xa <- rnorm(n*q)
dim(Xa) <- c(n,q)
Xb <- rnorm(n*p)
dim(Xb) <- c(n,p)


tildeUa <- MASS::Null(Xa)
XbXbT <- Xb %*% t(Xb)
XbStarXbStarT <- t(tildeUa) %*% XbXbT %*% tildeUa
tmp1 <- - tildeUa %*% solve(XbStarXbStarT) %*% t(tildeUa)
A <- tmp1 + mean(diag(tmp1))/(n-q) * tildeUa %*% t(tildeUa)

# generate reference normal random variables
offdiagA <- A - diag(diag(A))
diagNorm <- sqrt(sum(diag(A)^2))
ref1 <- NULL
ref2 <- NULL
for(i in 1:1000){
    ref1[i] <- diagNorm * rnorm(1)
    tmp <- rnorm(n)
    ref2[i] <- as.numeric(t(tmp) %*% offdiagA %*%tmp)
}


myResult <- NULL
jjj <- NULL
for(i in 1:10000){
    y <- rnorm(n)
    # proposed statistic
    theNumerator <- as.numeric( 
        t(y) %*% tmp1 %*% y
        )
    theDenominator <- sum(
                            (t(tildeUa) %*% y)^2
                          )
    theDenominator <- n-q
    proposedStat <- theNumerator / theDenominator
    
    
    # estimation of tau square    
    tildePa <- tildeUa %*% t(tildeUa)
    tildeEpsilon <- tildePa %*% y
    tauSquare <- 
        (
            (n-q)^2 * sum(tildeEpsilon^4)/sum(tildeEpsilon^2)^2-
                3*sum(diag(tildePa^2))
            )/
        (
            sum(tildePa^4)
            )+
        2
    tau <- sqrt(tauSquare)
    
    # referece distribution
    myRef <- tau * ref1 + ref2
    
    if(mean(myRef > ((n-q)*proposedStat)-sum(diag(tmp1))) <0.05){
        myResult[i] <- 1
    }
    else{
        myResult[i] <- 0
    }
}

mean(myResult)



