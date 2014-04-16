library(ks)
set.seed(1)
par(mfrow=c(2,1))

x<-rlnorm(100)
hist(x, col="red", freq=F)
lines(density(x))


y <- rkde(fhat=kde(x=x, h=hpi(x)), n=100, positive=TRUE)
hist(y, col="green", freq=F)

par(mfrow=c(1,1))
set.seed(1)
x<-rlnorm(100)
hist(x, prob=TRUE)

lines(density(x), col='red')

library(ks)
tmp <- kde(x, hpi(x))
lines(tmp$eval.points, tmp$estimate, col='green')

library(logspline)
lsfit <- logspline(x, lbound=0)
curve( dlogspline(x,lsfit), add=TRUE, col='blue' )

curve( dlnorm, add=TRUE, col='orange' )


