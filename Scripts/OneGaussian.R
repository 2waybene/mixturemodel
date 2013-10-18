setwd("/Users/li11/myGit/FRCompete/mixtureModels/output/")
# efg, 23 Aug 2006
# Study single Gaussian

pdf("SingleGaussian.pdf", width=8.0, height=10.0)

# Here P is the same as dnorm (probability density function for normal
# distribution), but other functions could be tried here.
P <- function(x, mean, sd)
{
  variance <- sd^2
  exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance)
}


# integrate(P, -1, 1, mean=0, sd=1) is same as integrate(dnorm, -1, 1, mean=0, sd=1)
# integrate(P,-5,5, mean=0, sd=1)   # should be close to 1.0

# Find "peaks" in array.
# R equivalent of Splus peaks() function.
# http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33097.html
# (see efg's posting to R-Help on 8 Feb 2007 about problem with ties.)
#
# peaks(c(1,4,4,1,6,1,5,1,1),3)
# [1] FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE

peaks <- function(series,span=3)
{
  z <- embed(series, span)
  s <- span%/%2
  v <- max.col(z, "first") == 1 + s   # take first if a tie
  result <- c(rep(FALSE,s),v)
  result <- result[1:(length(result)-s)]
  result
}


# First derivative.  Adjust x values to be center of interval.
# Spacing of x-points need not be uniform
Deriv1 <- function(x,y)
{
  y.prime <- diff(y) / diff(x)
  x.prime <- x[-length(x)] + diff(x)/2
  list(x = x.prime,
       y = y.prime)
}

# "Centered" 2nd-derivative. Spacing of x-points assumed to be uniform.
Deriv2 <- function(x,y)
{
  h <- x[2] - x[1]
  Range <- 2:(length(x)-1)  # Drop first and last points
  list(x = x[Range],
       y = (y[Range+1] - 2*y[Range] + y[Range-1]) / h^2)
}


# Setup test Gaussian with mean 0.0 and standard deviation 1.0
Delta <- 0.01
x <- seq(-5,5, by=Delta)
y <- P(x, 0.0, 1.0)

par(mfrow=c(3,1), oma=c(2,2,2,2))

##### Gaussian
plot(x,y, type="l", lwd=3,
     main="Single Gaussian",
     xlab="x", ylab="y")
abline(h=0, v=0, lty="dotted")

##### 1st Derivative
# According to Abramowitz & Stegun, inflection points are at +/- sigma.
derivative1 <- Deriv1(x,y)

plot(derivative1$x,derivative1$y, type="l",
     main="1st Derivative", xlab="x", ylab="y'")
abline(h=0, v=0, lty="dotted")

peaks.Deriv1   <- peaks( derivative1$y, span=3)
valleys.Deriv1 <- peaks(-derivative1$y, span=3)

points( derivative1$x[peaks.Deriv1], derivative1$y[peaks.Deriv1],
        pch=19, col="red")
points( derivative1$x[valleys.Deriv1], derivative1$y[valleys.Deriv1],
        pch=19, col="blue")

# Approximate location of peak and valley
derivative1$x[peaks.Deriv1]
derivative1$x[valleys.Deriv1]

s.approx <- (derivative1$x[valleys.Deriv1][1] - derivative1$x[peaks.Deriv1][1])/2
s.approx

#### 2nd Derivative
derivative2 <- Deriv2(x,y)

plot(derivative2$x,derivative2$y, type="l",
     main="2nd Derivative", xlab="x", ylab="y''")
abline(h=0, v=0, lty="dotted")

peaks.Deriv2   <- peaks( derivative2$y, span=3)
valleys.Deriv2 <- peaks(-derivative2$y, span=3)

points( derivative2$x[peaks.Deriv2], derivative2$y[peaks.Deriv2],
        pch=19, col="red")
points( derivative2$x[valleys.Deriv2], derivative2$y[valleys.Deriv2],
        pch=19, col="blue")

# Approximate location of mean of normal distribution:
derivative2$x[valleys.Deriv2]
derivative2$y[valleys.Deriv2]

mu.approx <-  derivative2$x[valleys.Deriv2]
mu.approx

# Peaks
derivative2$x[peaks.Deriv2]
derivative2$y[peaks.Deriv2]

mtext("U:/efg/lab/R/2006/MixturesOfDistributions/SingleGaussian.R",
      BOTTOM<-1, cex=0.8, adj=0.0, outer=TRUE)

# Fit distribution
# STUPID package:  cannot solve "exact case" without error.
# MUST have noise in data.
#fit1 <- nls(y ~ (a/b)*exp(-(x-c)^2/(2*b^2)),
#  start=list(a=1/sqrt(2*pi*s.approx^2), b=s.approx, c=mu.approx),
#  trace=TRUE)

# add small random noise (in general, how much noise is "enough"?)
set.seed(17)
y <- y + rnorm(length(y), 1E-3, 1E-4)
fit2 <- nls(y ~ (a/b)*exp(-(x-c)^2/(2*b^2)),
            start=list(a=1/sqrt(2*pi*s.approx^2), b=s.approx, c=mu.approx),
            control=nls.control(tol=1E-5, minFactor=1/1024),
            trace=TRUE)
fit2

1/sqrt(2*pi)


# Crude approximations
#derivative3 <- Deriv2(derivative1$x, derivative1$y)
#plot(derivative3$x,derivative3$y, type="l",
#  main="3rd Derivative", xlab="x", ylab="y''")
#abline(h=0, v=0, lty="dotted")
#
#derivative4 <- Deriv2(derivative2$x, derivative2$y)
#plot(derivative4$x,derivative4$y, type="l",
#  main="4th Derivative", xlab="x", ylab="y''")
#abline(h=0, v=0, lty="dotted")

dev.off()

#try a different Y now

y <- P(x, 1.1, 2.0)
set.seed(17)
y <- y + rnorm(length(y), 1E-3, 1E-4)
fit2 <- nls(y ~ (a/b)*exp(-(x-c)^2/(2*b^2)),
            start=list(a=1/sqrt(2*pi*s.approx^2), b=s.approx, c=mu.approx),
            control=nls.control(tol=1E-5, minFactor=1/1024),
            trace=TRUE)
fit2
