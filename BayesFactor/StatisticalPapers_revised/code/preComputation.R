# compute some quantities beforehand

# some quantities
tildeUa <- MASS::Null(Xa)
XbXbT <- Xb %*% t(Xb)
XbStarXbStarT <- t(tildeUa) %*% XbXbT %*% tildeUa
myBase <- - tildeUa %*% solve(XbStarXbStarT) %*% t(tildeUa)
tildePa <- tildeUa %*% t(tildeUa)
x0 <- sum(diag(myBase))/(n-q)
A <- myBase - x0 * tildePa

############################# GT Statistics
tmpX <- Xb
tmpX <- tmpX - Xa %*% solve(crossprod(Xa), crossprod(Xa, tmpX))
XX <- crossprod(t(tmpX))
var.num <- 2*sum(XX*XX)
lam <- eigen(XX, symmetric = TRUE, only.values=TRUE)$values
varGamma <- var(lam[1:(n-q)])

############################# EigenPrism: inference for high dimensional signal-to-noise ratios
ohoh <- eigen(Xb%*%t(Xb),symmetric = TRUE)
eigenXbXbT <- ohoh$values/p
myUb <- ohoh$vectors

myVariable <- Variable(n)
obj <- max(sum((myVariable^2) * eigenXbXbT^2), sum(myVariable^2))
myConstr <- list(sum(myVariable)==0, sum(eigenXbXbT * myVariable)==1)
myProblem <- Problem(Minimize(obj),myConstr)

ohMyTmp <- solve(myProblem)
myEPWeight <- ohMyTmp$getValue(myVariable)
# sometimes it returns "optimal_inaccurate"
#myValP1 <- ohMyTmp$value
# here is a little modification of the original method
#myValP1 <- max(sum(myEPWeight^2 * eigenXbXbT^2), sum(myEPWeight^2))
myValP1 <-  sum(myEPWeight^2)


######################### Xianyang Zhang
# deal with X
Gram <- t(X)%*%X/n
score.nodewiselasso <- getFromNamespace("score.nodewiselasso", "hdi")
node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                            parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
Theta <- node$out


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

