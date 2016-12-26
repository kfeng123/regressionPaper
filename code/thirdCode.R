
myExp=function(
    n=15,
    p=100,
    betaMod=0.00,
    varEpsilon=4,
    rho=runif(T,0,1),
    muX=runif(p,2,3),
    alpha=0,
    r=1,
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
            #     SigmaEst=var(t(X))
            #     myE=eigen(SigmaEst)$values
            #     jun=sum(myE)
            #     fang=sum(myE^2)-jun^2/(n-1)
            #     (temp-jun)/sqrt(fang)
        } 
        T=as.numeric(myTCal(X,y))
        
        
        ChenTCal=function(X,y){
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
            theTemp=0
            for(i1 in 1:n)for(i2 in 1:n){
                if(i1!=i2){
                    theTemp=theTemp+sum(X[,i1]*X[,i2])*y[i1]*y[i2]
                }
            }
            theTemp/ 
                (
                    t(y)%*%Q%*%y
                )
        }
        ChenT=as.numeric(ChenTCal(X,y))
        
        beta=0*beta
        Oh=NULL
        ChenOh=NULL
        for(ti in 1:100){
#            data=generateData(beta)
#            X=t(data$X)
#            y=data$y
            
#            TT=as.numeric(myTCal(X,y))
#            TT2=as.numeric(ChenTCal(X,y))
            
            myOrd=sample.int(n)
            TT=as.numeric(myTCal(X,y[myOrd]))
            TT2=as.numeric(ChenTCal(X,y[myOrd]))
            
            Oh=c(Oh,TT)
            ChenOh=c(ChenOh,TT2)
        }
        T=0+(mean(Oh>T)<=0.05) 
        ChenT=0+(mean(ChenOh>ChenT)<=0.05) 
        #        theTS=NULL
        #        for(th in 1:100) theTS[th]=myTCal(X,y[sample(n)])
        #        T=0+(mean(theTS>T)<=0.05)
        
        #        T=0+(T>qnorm(0.95))
        
        
        list(T=T,ChenT=ChenT,ChenTT=ChenTT)
    }
    
    REL=NULL
    ChenREL=NULL 
    ChenREL2=NULL
    for(i in 1:100){
        temp=simul()
        REL[i]=temp$T
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

myExp(n=10,p=100,r=1,TheLambda = 100,varEpsilon = 4,betaMod =0.01)


library(xtable)
myTable=xtable(resul)
digits(myTable)=c(0,0,0,2,2,2)
print(myTable,include.rownames = FALSE)