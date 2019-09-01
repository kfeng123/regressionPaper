# deal with X
Gram <- t(X)%*%X/n
score.nodewiselasso <- getFromNamespace("score.nodewiselasso", "hdi")
node <- score.nodewiselasso(X, wantTheta=TRUE, verbose=FALSE, lambdaseq="quantile",
                            parallel=FALSE, ncores=2, oldschool = FALSE, lambdatuningfactor = 1)
Theta <- node$out