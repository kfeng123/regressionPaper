
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
    TheLambda=p){
    
    beta=rep(c(1,-1),p/2)
    #beta=c(rep(1,p/2),rep(0,p/2))
    #beta=c(rep(1,5),rep(0,p-5))
    
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
        V=rnorm(p*r,0,1)
        dim(V)=c(p,r)
        V=svd(V)$u
        theSigma=rep(0,p*p)
        dim(theSigma)=c(p,p)
        theSigma=V%*%D%*%D%*%t(V)+diag(rep(1,p))
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
                 (temp-p*lambdaEst)/lambdaEst/sqrt(2*p)
        } 
        #T=0+(pnorm(myTCal(X,y))>=(1-0.05))
        T=myTCal(X,y)
        
        
        
        list(T=T)
    }
    
    REL=NULL
    for(i in 1:500){
        temp=simul()
        REL[i]=temp$T
    }
    
    #xxx=NULL
    #for(j in 1:length(REL)){
    #xxx[j]=qchisq((j-0.5)/length(REL),df=1)
    #}
    #plot(xxx,sort(REL))
    #abline(0,1)
    
    #hist(pchisq(REL,df=1))
    
    list(myPower=mean(pnorm(REL)>=1-0.05))
    
}

myExp(n=60,p=300,r=0,TheLambda = 100,varEpsilon = 1,betaMod =100)


library(xtable)
myTable=xtable(resul)
digits(myTable)=c(0,0,0,2,2,2)
print(myTable,include.rownames = FALSE)