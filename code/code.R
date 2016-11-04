n=80
p=550
T=10
betaMod=0.02
varEpsilon=4
rho=runif(T,0,1)
muX=runif(p,2,3)
alpha=1
beta=c(rep(1,p/2),rep(0,p/2))
#beta=c(rep(1,5),rep(0,p-5))
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
    list(X=X,y=y)
}

Q=diag(n)-rep(1,n)%*%t(rep(1,n))/n

simul=function(){
    data=generateData()
    X=t(data$X)
    y=data$y
    
    fenzi=(t(rep(1,n))%*%solve(t(X)%*%X)%*%Q%*%y)^2
    myVar=function(y){
        1/(n-2)*t(y)%*%Q%*%(
            diag(n)-
                solve(t(X)%*%X)%*%rep(1,n)%*%t(rep(1,n))%*%solve(t(X)%*%X)/
                as.numeric(t(rep(1,n))%*%solve(t(X)%*%X)%*%Q%*%solve(t(X)%*%X)%*%rep(1,n))
            )%*%Q%*%y
    }
    fenmu=myVar(y)*(t(rep(1,n))%*%solve(t(X)%*%X)%*%Q%*%solve(t(X)%*%X)%*%rep(1,n))
    
    T=fenzi/fenmu
}

REL=NULL

for(i in 1:1000){
    REL[i]=simul()
}

#xxx=NULL
#for(j in 1:length(REL)){
    #xxx[j]=qchisq((j-0.5)/length(REL),df=1)
#}
#plot(xxx,sort(REL))
#abline(0,1)

#hist(pchisq(REL,df=1))

(power=mean(pchisq(REL,df=1,lower.tail = FALSE)<0.05))

