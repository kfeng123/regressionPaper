mySeed=1
myExp=function(
    n=15,
    p=100,
    betaMod=0.00,
    varEpsilon=4,
    rho=runif(T,0,1),
    muX=runif(p,2,3),
    alpha=0,
    r=1,
    k=1,
    Alt=TRUE,
    TheLambda=p){
    
    #beta=rep(c(1,-1),p/2)
    #beta=c(rep(1,p/2),rep(0,p/2))
    beta=c(rep(1,5),rep(0,p-5))
    
    #beta=rnorm(p)
    #beta=beta[sample.int(p)]
    
    beta=beta/sqrt(sum(beta^2))
    beta=sqrt(betaMod)*beta
    
    if(r==0){
        theSigma=diag(rep(1,p))
    }else{
        if(r!=1){
            D=diag(rep(sqrt(TheLambda),r))
        }else{
            D=TheLambda
            dim(D)=c(1,1)
        }
        set.seed(mySeed)
        V=rnorm(p*r,0,1)
        dim(V)=c(p,r)
        V=svd(V)$u
        theSigma=rep(0,p*p)
        dim(theSigma)=c(p,p)
        theSigma=V%*%D%*%D%*%t(V)+diag(rep(1,p))
        if(Alt==FALSE){
            beta=(diag(p)-V%*%t(V))%*%beta
            beta=beta/sqrt(sum(beta^2))
            beta=sqrt(betaMod)*beta
        }
    }
    
    
    #    beta=(diag(p)-V%*%t(V))%*%rnorm(p)
    #    beta=beta/sqrt(sum(beta^2))
    #    beta=sqrt(betaMod)*beta
    
    
    
    generateData=function(beta){
        X=rep(0,n*p)
        dim(X)=c(n,p)
        if(r==0){
            Z=rnorm(n*p,0,1)
            dim(Z)=c(n,p)
            X=Z
        }else{
            U=rnorm(n*r,0,1)
            dim(U)=c(n,r)
            Z=rnorm(n*p,0,1)
            dim(Z)=c(n,p)
            X=U%*%D%*%t(V)+Z
        }
        y=X%*%beta+rep(alpha,n)+rnorm(n,0,sqrt(varEpsilon))
        list(X=X,y=y)
    }
    
    Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
    W=eigen(Q)$vectors[,1:(n-1)]
    
    
    ChenTT=pnorm(qnorm(0.05)+n*sum((theSigma%*%beta)^2)/sqrt(2*sum(theSigma^2))/4)
    
    
    simul=function(){
        data=generateData(beta)
        X=t(data$X)
        y=data$y
        
        myTCal=function(X,y){
            temp=(
                sum(y^2)-
                    n*(mean(y))^2
            )/ 
                (
                    t(y)%*%W%*%solve(t(W)%*%t(X)%*%X%*%W)%*%t(W)%*%y
                )
                 myE=eigen(t(W)%*%t(X)%*%X%*%W)$values
                 lambdaEst=sum(myE[(k+1):(n-1)])/(p-k)/(n-1)
                 lambdaEst=1
                 (temp-p*lambdaEst)/lambdaEst/sqrt(2*p)
        } 
        #T=0+(pnorm(myTCal(X,y))>=(1-0.05))
        T=myTCal(X,y)
        
        pcaX=prcomp(t(X))$x[,k]
        dim(pcaX)=c(n,k)
        U=cbind(1,pcaX)
        pE=solve(t(U)%*%U)%*%t(U)%*%y
        pE2=pE[-1]
        dim(pE2)=c(k,1)
        fenzi=t(pE2)%*%solve(solve(t(U)%*%U)[2:(k+1),2:(k+1)])%*%pE2/k
        fenmu=t(y)%*%(diag(n)-U%*%solve(t(U)%*%U)%*%t(U))%*%y/(n-k-1)
        Fstat=fenzi/fenmu
        
        list(T=T,Fstat=Fstat,ChenTT=ChenTT)
    }
    
    REL=NULL
    FEL=NULL
    CEL=NULL
    for(i in 1:500){
        temp=simul()
        REL[i]=temp$T
        FEL[i]=temp$Fstat
        CEL[i]=temp$ChenTT
    }
    
    
    list(myPower=mean(pnorm(REL)>=(1-0.05)),
         FPower=mean(FEL>=qf(1-0.05,k,n-k-1)),
    CPower=mean(ChenTT))
    
}

myExp(n=150,p=500,r=1,k=1, Alt=TRUE,TheLambda = sqrt(500),varEpsilon = 1,betaMod =1)


library(xtable)
myTable=xtable(resul)
digits(myTable)=c(0,0,0,2,2,2)
print(myTable,include.rownames = FALSE)