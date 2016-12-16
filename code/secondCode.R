
myExp=function(
    n=20,
    p=100,
    T=1,
    betaMod=0.00,
    varEpsilon=4,
    rho=runif(T,0,1),
    muX=runif(p,2,3),
    alpha=0){
    
    #beta=rep(c(1,-1),p/2)
    beta=c(rep(1,p/2),rep(0,p/2))
    #beta=c(rep(1,5),rep(0,p-5))
    
    beta=beta[sample.int(p)]
     
    beta=beta/sqrt(sum(beta^2))
    beta=sqrt(betaMod)*beta
    
    
    theSigma=rep(0,p*p)
    dim(theSigma)=c(p,p)
    for(i in 1:p)for(j in 1:p){
        if(abs(i-j)<T){
            temp=0
            for(k in 1:(T-abs(i-j))){
                temp=temp+rho[k]*rho[k+abs(i-j)]
            }
            theSigma[i,j]=temp
        }
    }
    theSigma=theSigma+10000
    
    generateData=function(beta){
        X=rep(0,n*p)
        dim(X)=c(n,p)
        Z=rnorm(n*(p+T-1),0,1)
        dim(Z)=c(n,p+T-1)
        
        for(i in 1:n)for(j in 1:p){
            X[i,j]=muX[j]+sum(rho*Z[i,j:(j+T-1)])
        }
        
        iit=rep(rnorm(n,0,100),each=p)
        dim(iit)=c(p,n)
        
        X=X+t(iit)
        y=X%*%beta+rep(alpha,n)+rnorm(n,0,sqrt(varEpsilon))
        
        list(X=X,y=y)
    }
    
    Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
    W=eigen(Q)$vectors[,1:(n-1)]
    
    
    ChenT=pnorm(qnorm(0.05)+n*sum((theSigma%*%beta)^2)/sqrt(2*sum(theSigma^2))/4)
    
    
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
        
        beta=0*beta
        Oh=NULL
        for(ti in 1:50){
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
            TT=as.numeric(myTCal(X,y))
            Oh=c(Oh,TT)
        }
        T=0+(mean(Oh>T)<=0.05) 
#        theTS=NULL
#        for(th in 1:100) theTS[th]=myTCal(X,y[sample(n)])
#        T=0+(mean(theTS>T)<=0.05)
        
#        T=0+(T>qnorm(0.95))
        
        #thePhi=function(i1,i2,i3,i4){
        #    1/4*t(X[,i1]-X[,i2])%*%(X[,i3]-X[,i4])*(y[i1]-y[i2])*(y[i3]-y[i4])
        #}
        #theTemp=0
        #for(i1 in 1:n)for(i2 in 1:n)for(i3 in 1:n)for(i4 in 1:n){
        #    if(i1!=i2&i1!=i3&i1!=i4&i2!=i3&i2!=i4&i3!=i4){
        #        theTemp=theTemp+thePhi(i1,i2,i3,i4)
        #    }
        #}
        #theTemp=theTemp*n*(n-1)*(n-2)*(n-3)/4/3/2/1
        #ChenT=n*theTemp/sqrt(2*sum(theSigma^2))/4
        
        list(T=T,ChenT=ChenT)
    }
    
    REL=NULL
    ChenREL=NULL 
    for(i in 1:100){
        temp=simul()
        REL[i]=temp$T
        ChenREL[i]=temp$ChenT
    }
    
    #xxx=NULL
    #for(j in 1:length(REL)){
    #xxx[j]=qchisq((j-0.5)/length(REL),df=1)
    #}
    #plot(xxx,sort(REL))
    #abline(0,1)
    
    #hist(pchisq(REL,df=1))
    
    list(myPower=mean(REL),
    chenPower=mean(ChenREL))
    
}

myExp(varEpsilon = 4,betaMod =0.000001,T=1 )

jjj=NULL
 for(myMod in c(0,0.04)){
     jjj=c(jjj,list(list(n=40,p=310,myMod=myMod)))
 }
 for(myMod in c(0,0.04)){
     jjj=c(jjj,list(list(n=80,p=550,myMod=myMod)))
 }
 
 outN=NULL
 outP=NULL
 outMod=NULL
 outChen=NULL
 outNew=NULL
 for(pa in jjj){
     outN=c(outN,pa$n)
     outP=c(outP,pa$p)
     outMod=c(outMod,pa$myMod)
     tmp=myExp(n=pa$n,p=pa$p,betaMod=pa$myMod)
     outChen=c(outChen,tmp$chenPower)
     outNew=c(outNew,tmp$myPower)
 }
 resul=data.frame(outN,outP,outMod,outChen,outNew)
 
 library(xtable)
 myTable=xtable(resul)
 digits(myTable)=c(0,0,0,2,2,2)
 print(myTable,include.rownames = FALSE)