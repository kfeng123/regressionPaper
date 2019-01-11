library(MASS)
source("./gt.R")

# RepTime <- 2000
# M <- 2000
# alpha <- 0.05
# betabGen <- "dense"
# XGen <- "equalCor"
# epsilonDis <- "t"
# 
# n <- 100
# q <- 10
# p <- 1000

# data generation
Xa <- rnorm(n*q)
dim(Xa) <- c(n,q)
if(XGen == "iidnormal"){
    Xb <- rnorm(n*p)
    dim(Xb) <- c(n,p)
}
if(XGen == "equalCor" ){
    Xb <- sqrt(0.9)*rnorm(n*p)
    dim(Xb) <- c(n,p)
    for(i in 1:n){
        Xb[i,] <- Xb[i,] + sqrt(0.1)*rnorm(1)
    }
}
if(XGen == "Toeplitz"){
    tmp <- matrix(rep(0,p*p),p)
    for(i in 1:p) for (j in 1:p)
        tmp[i,j] <- 0.9^(abs(i-j))
    Xb <- mvrnorm(n,rep(0,p),tmp)
}
if(XGen == "real"){
    Data <- read.csv("riboflavin.csv")
    rownames(Data)
    Data<- Data[,-1]
    Data <- t(as.matrix(Data))
    X <- Data[,-1]
    y <- Data[,1]
    Xa <- X[,1:10]
    Xb <- X[,-(1:10)]
    n <- nrow(Xa)
    q <- ncol(Xa)
    p <- ncol(Xb)
}




youX <- cbind(Xa,Xb)
colnames(youX) <- seq(1,p+q)


### GT Statistics

tmpX <- Xb
tmpX <- tmpX - Xa %*% solve(crossprod(Xa), crossprod(Xa, tmpX))
XX <- crossprod(t(tmpX))
var.num <- 2*sum(XX*XX)
lam <- eigen(XX, symmetric = TRUE, only.values=TRUE)$values
varGamma <- var(lam[1:(n-q)])



# some quantities
tildeUa <- MASS::Null(Xa)
XbXbT <- Xb %*% t(Xb)
XbStarXbStarT <- t(tildeUa) %*% XbXbT %*% tildeUa
myBase <- - tildeUa %*% solve(XbStarXbStarT) %*% t(tildeUa)
tildeUatildeUaT <- tildeUa %*% t(tildeUa) 
x0 <- sum(diag(myBase))/(n-q)
A <- myBase - x0 * tildeUatildeUaT
tildePa <- tildeUa %*% t(tildeUa)

# generate reference normal random variables
offdiagA <- A - diag(diag(A))
diagNorm <- sqrt(sum(diag(A)^2))
ref1 <- NULL
ref2 <- NULL
for(i in 1:M){
    ref1[i] <- diagNorm * rnorm(1)
    tmp <- rnorm(n)
    ref2[i] <- as.numeric(t(tmp) %*% offdiagA %*%tmp)
}

# # Base
# offdiagBase <- myBase - diag(diag(myBase))
# diagBase <- diag(myBase)
# refBase1 <- NULL
# refBase2 <- NULL
# # Plus
# offdiagPlus <- tildeUatildeUaT - diag(diag(tildeUatildeUaT))
# diagPlus <- diag(tildeUatildeUaT)
# refPlus1 <- NULL
# refPlus2 <- NULL
# for(i in 1:M){
#     tmp <- rnorm(n)
#     refBase1[i] <- sum(diagBase * tmp)
#     refPlus1[i] <- sum(diagPlus * tmp)
#     tmp <- rnorm(n)
#     refBase2[i] <- as.numeric(t(tmp) %*% offdiagBase %*% tmp)
#     refPlus2[i] <- as.numeric(t(tmp) %*% offdiagPlus %*% tmp)
# }

outMy <- NULL
outGt <- NULL
for(SNR in c(0,5,10,15,20,25,30)){
    myResultOld <- NULL
    gtResult <- NULL
    for(i in 1:RepTime){
        
        # betaGen
        if(betabGen == "dense") {
            betabO <- runif(p,-1,1)
        }
        if(betabGen == "sparse") {
            betabO <- runif(p,-1,1)
            betabO[sample(p,p/20*19)] <- 0
        }
            
        meanSig <- Xb%*% betabO
            
        # Generate y
        tmpSNR <- sqrt((n-q)*varGamma)*sum(betabO^2)/p
        #betab <- betabO/sqrt(tmpSNR)*sqrt(SNR)
        if(epsilonDis == "t")
            innov <- rt(n,8)
        if(epsilonDis == "chi")
            innov <- (rchisq(n,4)-4)/sqrt(8)
        y <- innov + meanSig/sqrt(tmpSNR)*sqrt(SNR)
        
        
        ### GT statistics
        Y <- y
        S <- sum(Y * (XX %*% Y)) / sum((t(tildeUa) %*% Y)^2)
        
        lams <- lam
        lams[1:(n-q)] <- lams[1:(n-q)] - S
        
        p.value <- .getP(lams)
        
        gtResult[i] <- (p.value <alpha)
    
                                
        # proposed statistic
        theNumerator <- as.numeric( 
            t(y) %*% myBase %*% y
            )
        theDenominator <- sum(
                                (t(tildeUa) %*% y)^2
                              )
        proposedStat <- theNumerator / theDenominator
        
        
        # estimation of tau square    
        tildeEpsilon <- tildePa %*% y
        tauSquare <- 
            (
                (n-q)^2 * sum(tildeEpsilon^4)/sum(tildeEpsilon^2)^2-
                    3*sum(diag(tildePa)^2)
                )/
            (
                sum(tildePa^4)
                )+
            2
        if(tauSquare < 0) tauSquare <- 0
        tau <- sqrt(tauSquare)
        
        # reference distribution
        myRef <- tau * ref1 + ref2
        if(mean(myRef > ((n-q)*proposedStat)-sum(diag(myBase))) <alpha){
            myResultOld[i] <- 1
        }
        else{
            myResultOld[i] <- 0
        }
        
        # myRefStep1 <- tau * (refBase1 - x0 * refPlus1) + refBase2 - x0 * refPlus2
        # x1 <- quantile(myRefStep1,1-alpha)/(n-q) + x0
        # 
        # myRefStep2 <- tau * (refBase1 - x1 * refPlus1) + refBase2 - x1 * refPlus2
        # x2 <- quantile(myRefStep2,1-alpha)/(n-q) + x0
        # if(proposedStat > x2){
        #     myResult[i] <- 1
        # }
        # else{
        #     myResult[i] <- 0
        # }
        # 
        # 
        # if(proposedStat > x1){
        #     myResultTmp[i] <- 1
        # }
        # else{
        #     myResultTmp[i] <- 0
        # }
    }
    
    outMy <- c(outMy,mean(myResultOld))
    outGt <- c(outGt,mean(gtResult))
}
ohResult <- data.frame('SNR'=c(0,5,10,15,20,25,30),'outMy'=outMy,'outGt'=outGt)



