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

if(XGen == "uniform"){
    tmp <- rnorm(p*p)
    dim(tmp) <- c(p,p)
    tmpp <- tmp%*% t(tmp)
    Xb <- mvrnorm(n,rep(0,p), tmpp)
}

if(XGen == "Toeplitz"){
    tmp <- matrix(rep(0,p*p),p)
    for(i in 1:p) for (j in 1:p)
        tmp[i,j] <- 0.9^(abs(i-j))
    Xb <- mvrnorm(n,rep(0,p),tmp)
}

if(XGen == "factor"){
    K <- 2
    factorB <- rnorm(p*K)
    dim(factorB) <- c(p,K)
    factorF <- rnorm(n*K)
    dim(factorF) <- c(n,K)
    factorU <- rnorm(n*p)
    dim(factorU) <- c(n,p)
    
    
    factorU <- sqrt(0.9)*rnorm(n*p)
    dim(factorU) <- c(n,p)
    for(i in 1:n){
        factorU[i,] <- factorU[i,] + sqrt(0.1)*rnorm(1)
    }
    
    
    Xb <- factorF%*%t(factorB) + factorU
    
    tmpX <- factorU
    tmpX <- tmpX - Xa %*% solve(crossprod(Xa), crossprod(Xa, tmpX))
    XX <- crossprod(t(tmpX))
    var.num <- 2*sum(XX*XX)
    lam <- eigen(XX, symmetric = TRUE, only.values=TRUE)$values
    factorVarU <- var(lam[1:(n-q)])
}

if(XGen == "real"){
    Data <- read.csv("riboflavin.csv")
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

if(XGen == "rowCor"){
    Xb <- rnorm(n*p)
    dim(Xb) <- c(n,p)
    tmpMatrix <- rep(0.1,n*n)
    dim(tmpMatrix) <- c(n,n)
    diag(tmpMatrix) <- 1
    tmpChol <- Matrix::chol(tmpMatrix)
    Xb <- t(tmpChol) %*% Xb
}


youX <- cbind(Xa,Xb)
colnames(youX) <- seq(1,p+q)