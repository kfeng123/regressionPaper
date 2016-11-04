myExp=function(
    n=40,
    p=310,
    T=20,
    betaMod=0.06,
    varEpsilon=4,
    rho=runif(T,0,1),
    muX=runif(p,2,3),
    alpha=1){
    
    #beta=c(rep(1,p/2),rep(0,p/2))
    beta=c(rep(1,5),rep(0,p-5))
    
    beta=beta[sample.int(p)]
     
    beta=beta/sqrt(sum(beta^2))
    beta=sqrt(betaMod)*beta
    generateData=function(){
        X=rep(0,n*p)
        dim(X)=c(n,p)
        Z=rnorm(n*(p+T-1),0,1)
        dim(Z)=c(n,p+T-1)
        
        for(i in 1:n)for(j in 1:p){
            X[i,j]=muX[j]+sum(rho*Z[i,j:(j+T-1)])
        }
        y=X%*%beta+rep(alpha,n)+rnorm(n,0,sqrt(varEpsilon))
        
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
        list(X=X,y=y,theSigma=theSigma)
    }
    
    Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n
    
    simul=function(){
        data=generateData()
        X=t(data$X)
        y=data$y
        theSigma=data$theSigma
        
        myInv=solve(t(X)%*%X)
        fenzi=(t(rep(1,n))%*%myInv%*%Q%*%y)^2
        myVar=function(y){
            1/(n-2)*t(y)%*%Q%*%(
                diag(n)-
                    solve(t(X)%*%X)%*%rep(1,n)%*%t(rep(1,n))%*%solve(t(X)%*%X)/
                    as.numeric(t(rep(1,n))%*%solve(t(X)%*%X)%*%Q%*%solve(t(X)%*%X)%*%rep(1,n))
                )%*%Q%*%y
        }
        fenmu=myVar(y)*(t(rep(1,n))%*%myInv%*%Q%*%myInv%*%rep(1,n))
        
        T=fenzi/fenmu
        
        ChenT=pnorm(qnorm(0.05)+n*sum((theSigma%*%beta)^2)/sqrt(2*sum(theSigma^2))/4)
        
        
        
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
    for(i in 1:200){
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
    
    list(myPower=mean(pchisq(REL,df=1,lower.tail = FALSE)<0.05),
    #chenPower=mean(pnorm(ChenREL,lower.tail = FALSE)<0.05))
    chenPower=mean(ChenREL))
    
}


jjj=NULL
for(myMod in c(0,0.02,0.04,0.06)){
    jjj=c(jjj,list(list(n=40,p=310,myMod=myMod)))
}
for(myMod in c(0,0.02,0.04,0.06)){
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
