lassoTest.statistic <- max(abs(sqrt(n)*betabO_real + lassoTest.leftMultilier %*% innov / sqrt(n) )) 

if( mean( lassoTest.reference > lassoTest.statistic ) < alpha ){
    lassoTest_one_result <- 1
}else{
    lassoTest_one_result <- 0
}
