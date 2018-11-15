myTCal=function(X,y){
    n <- nrow(X)
    p <- ncol(X)
    kappa <- sum(X^2)/n/p
    #temp <- (
            #t(y)%*%X%*%solve(t(X)%*%X+diag(rep(kappa,p))) %*% t(X) %*%y
    #)/ 
        #(
            #sum(y^2)
        #)
    
    
    
    temp <- (
            sum(y^2)
    )/ 
        (
            t(y)%*%solve (X %*% t(X)) %*%y
        )
} 



ChenTCal=function(X,y){
    n <- nrow(X)
    p <- ncol(X)
    # thePhi=function(i1,i2,i3,i4){
    #     1/4*t(X[,i1]-X[,i2])%*%(X[,i3]-X[,i4])*(y[i1]-y[i2])*(y[i3]-y[i4])
    # }
    # theTemp=0
    # for(i1 in 1:n)for(i2 in 1:n)for(i3 in 1:n)for(i4 in 1:n){
    #     if(i1!=i2&i1!=i3&i1!=i4&i2!=i3&i2!=i4&i3!=i4){
    #         theTemp=theTemp+thePhi(i1,i2,i3,i4)
    #     }
    # }
    # theTemp=theTemp*n*(n-1)*(n-2)*(n-3)/4/3/2/1
    # ChenT=n*theTemp/sqrt(2*sum(theSigma^2))/4
    theTemp <- 0
    for(i1 in 1:n)for(i2 in 1:n){
        if(i1!=i2){
            theTemp=theTemp+sum(X[i1,]*X[i2,])*y[i1]*y[i2]
        }
    }
    theTemp/ 
        (
            sum(y^2)
        )
}



myExp <- function(n=15, p=100, betaMod=1){
    
    #beta=c(rep(1,p/2),rep(0,p/2))
    #beta=c(rep(1,5),rep(0,p-5))
    
    #beta=beta[sample.int(p)]
    
    
    X <- rnorm(n*p,0,1)
    dim(X) <- c(n,p)
    generateData=function(beta){
        
        y=X%*%beta+rnorm(n,0,1)
        
        list(X=X,y=y)
    }
    
    ChenTT <- pnorm(qnorm(0.05)+n*betaMod/sqrt(2*p))

    simul <- function(){
        #beta=rep(c(1,-1),p/2)
        beta=rnorm(p)
        beta=beta/sqrt(sum(beta^2))
        beta=sqrt(betaMod)*beta
        data <- generateData(beta)
        X <- data$X
        y <- data$y
        
        # Bayes
        kappa <- sum(X^2)/n/p
        haha <- X%*%solve(t(X)%*%X+diag(rep(kappa,p))) %*% t(X)
        temp <- (
                t(y)%*% haha %*%y
        )/ 
            (
                sum(y^2)
            )
        
        myT <- as.numeric(temp)
        
        ChenT <- as.numeric(ChenTCal(X,y))
        
        Oh <- rep(0,100)
        ChenOh <- rep(0,100)
        for(ti in 1:100){
            myOrd=sample.int(n)
            temp <- (
                    t(y[myOrd])%*% haha %*%y[myOrd]
            )/ 
                (
                    sum(y^2)
                )
            TT=as.numeric(temp)
            TT2=as.numeric(ChenTCal(X,y[myOrd]))
            
            Oh[ti] <- TT
            ChenOh[ti] <- TT2
        }
        myT <- 0+(mean(Oh>myT)<=0.05) 
        ChenT <- 0+(mean(ChenOh>ChenT)<=0.05) 
        #        theTS=NULL
        #        for(th in 1:100) theTS[th]=myTCal(X,y[sample(n)])
        #        T=0+(mean(theTS>T)<=0.05)
        
        #        T=0+(T>qnorm(0.95))
        
        
        list(myT=myT,ChenT=ChenT,ChenTT=ChenTT)
    }
    
    REL=NULL
    ChenREL=NULL 
    ChenREL2=NULL
    for(i in 1:300){
        temp=simul()
        REL[i]=temp$myT
        ChenREL[i]=temp$ChenT
        ChenREL2[i]=temp$ChenTT
    }
    
    #xxx=NULL
    #for(j in 1:length(REL)){
    #xxx[j]=qchisq((j-0.5)/length(REL),df=1)
    #}
    #plot(xxx,sort(REL))
    #abline(0,1)
    
    #hist(pchisq(REL,df=1))
    
    list(myPower=mean(REL),
         chenPower=mean(ChenREL),ChenTT=mean(ChenREL2))
    
}

(resul <- myExp(n=30,p=60,betaMod =5))


#library(xtable)
#myTable=xtable(resul)
#digits(myTable)=c(0,0,0,2,2,2)
#print(myTable,include.rownames = FALSE)