 library(BigTSP)

 library(randomForest)
 x=matrix(rnorm(100*20),100,20)
 y=rbinom(100,1,0.5)
 y=as.factor(y)
 fit=tsp.randomForest(x,y)
 predict(fit,x[1:10,])

 plot(fit, main="A simulated case")
