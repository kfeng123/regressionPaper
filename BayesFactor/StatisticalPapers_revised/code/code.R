#RepTime <- 2000
#M <- 2000
#alpha <- 0.05
#betabGen <- "dense"
#XGen <- "equalCor"
#epsilonDis <- "t"
#
#n <- 100
#q <- 10
#p <- 1000

source("./gt.R")
SNRsequence <- c(0,10,20,30)
source("./dataGeneration.R")
source("./preComputation.R")
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
outEp <- NULL
outSNR <- NULL
for(SNR in SNRsequence){
    myPb <- txtProgressBar(style = 3)
    for( j in 1:plotPointSize ){
        myResultOld <- NULL
        gtResult <- NULL
        epResult <- NULL
            
        # betaGen
        if(betabGen == "dense") {
            if(betaDistribution == "unif")
                betabO <- runif(p,-1,1)
            if(betaDistribution == "normal")
                betabO <- rnorm(p)
        }
        if(betabGen == "sparse") {
            if(betaDistribution == "unif")
                betabO <- runif(p,-1,1)
            if(betaDistribution == "normal")
                betabO <- rnorm(p)
            betabO[sample(p,p/20*19)] <- 0
        }
            
        meanSig <- Xb%*% betabO
        if(XGen == "factor"){
            meanSig <- factorU %*% betabO
        }
        
        for(i in 1:RepTime){
                
            # Generate y
            #betab <- betabO/sqrt(tmpSNR)*sqrt(SNR)
            if(epsilonDis == "t"){
                myTdf <- 9
                innov <- rt(n,myTdf)
                myPhi <- (myTdf-2)/myTdf
            }
            if(epsilonDis == "chi"){
                innov <- (rchisq(n,4)-4)/sqrt(8)
                myPhi <- 1
            }
            tmpSNR <- sqrt((n-q)*varGamma)*myPhi*sum(betabO^2)/p
            if(XGen == "factor"){
                tmpSNR <- sqrt((n-q)*factorVarU)*myPhi*sum(betabO^2)/p
            }
            y <- innov + meanSig/sqrt(tmpSNR)*sqrt(SNR)
            
            
            # EigenPrism: inference for high dimensional signal-to-noise ratios
            myZ <- as.numeric( t(myUb) %*% y )
            ohT <- sum(myEPWeight * myZ^2 )
            epResult[i] <- ( ohT/ ( sqrt(2*myValP1) * sum(y^2)/n ) > qnorm(1-alpha) )
            
            
            
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
        outEp <- c(outEp,mean(epResult))
        outSNR <- c(outSNR, SNR)
        setTxtProgressBar(myPb,j/plotPointSize)
    }
    close(myPb)
}
ohResult <- data.frame('SNR'=outSNR,'outMy'=outMy,'outGt'=outGt, 'outEp' = outEp)

write.csv(ohResult,paste0("tmp/",n,"_",XGen,"_",epsilonDis,"_",betabGen,".csv"),row.names=FALSE)



