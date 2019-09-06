# X
# y
# set
LT.M <- 500
here_alpha <- 0.95
set <- 1:q

sreg <- scalreg(X, y)
beta.hat <- sreg$coefficients
sigma.sq <- sum((y-X%*%beta.hat)^2)/(n-sum(abs(beta.hat)>1e-9))
beta.db <- beta.hat+Theta%*%t(X)%*%(y-X%*%beta.hat)/n
Omega <- diag(Theta%*%Gram%*%t(Theta))*sigma.sq
stat.boot.st <- stat.boot.nst <- rep(NA,LT.M)

for (xianyang.i in 1:LT.M) {
    e <- rnorm(n)
    xi.boot <- Theta[set,]%*%t(X)%*%e*sqrt(sigma.sq)/sqrt(n)
    stat.boot.nst[xianyang.i] <- max(abs(xi.boot))
    stat.boot.st[xianyang.i] <- max(abs(xi.boot)/sqrt(Omega[set]))
}

crit.nst <- quantile(stat.boot.nst, here_alpha)
crit.st <- quantile(stat.boot.st, here_alpha)

up.nst <- beta.db[set] + crit.nst/sqrt(n)
low.nst <- beta.db[set] - crit.nst/sqrt(n)

band.nst <- rbind(low.nst, up.nst)

# test
if(prod(up.nst > 0) * prod(low.nst < 0)){
    # accept
    lassoTest_one_result <- 0
}else{
    # reject
    lassoTest_one_result <- 1
}
