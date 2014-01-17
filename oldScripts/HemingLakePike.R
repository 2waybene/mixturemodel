# efg, 15 Aug 2006
# Data from http://www.math.mcmaster.ca/peter/mix/demex/expike.html

#pdf("HemingLakePike.pdf", width=8.0, height=10.0)

library(Hmisc)

weighted.stats <- function(name, weights)
{
  mean <- wtd.mean(d$MID.POINT, weights)
  sd   <- sqrt(wtd.var(d$MID.POINT, weights))
  N    <- sum(weights)
  return( list(mean=mean,sd=sd, N=N) )
}

plot.normal.curve <- function(x, y, line.color)
{
  stats <- weighted.stats("1", y)
  y.norm <- sum(y) * mean(diff(d$MID.POINT)) * dnorm(x, stats$mean, stats$sd)
  lines(x, y.norm,
        lty="dotted", col=line.color)
  cat(stats$mean, stats$sd, stats$N, "\n")
  y.norm
}

d <- read.table("http://research.stowers-institute.org/efg/R/Statistics/MixturesOfDistributions/HemingLakePike.dat",
                header=TRUE)

par(mfrow=c(6,1), mar=c(2,2,2,2))
plot(d$MID.POINT, d$TOTAL, type="l", col="black", ylim=c(0,20),
     xaxt="none", main="Heming Lake Pike", lwd=2)

plot(d$MID.POINT, d$AGE_1, type="l", col="red", ylim=c(0,20),
     xaxt="none", main="1")
y1 <- plot.normal.curve(d$MID.POINT, d$AGE_1, "red")

plot(d$MID.POINT, d$AGE_2, type="l", col="green",ylim=c(0,20),
     xaxt="none", main="2")
y2 <- plot.normal.curve(d$MID.POINT, d$AGE_2, "green")

plot(d$MID.POINT, d$AGE_3, type="l", col="blue", ylim=c(0,20),
     xaxt="none", main="3")
y3 <- plot.normal.curve(d$MID.POINT, d$AGE_3, "blue")

plot(d$MID.POINT, d$AGE_4, type="l", col="magenta",ylim=c(0,20),
     xaxt="none", main="4")
y4 <- plot.normal.curve(d$MID.POINT, d$AGE_4, "magenta")

plot(d$MID.POINT, d$AGE_5, type="l", col="cyan", ylim=c(0,20),
     xaxt="none", main="5")
y5 <- plot.normal.curve(d$MID.POINT, d$AGE_5, "cyan")


# efg, 22 Aug 2006
# Data adapted from http://www.math.mcmaster.ca/peter/mix/demex/expike.html

# integrate(P, -1, 1, mean=0, sd=1) = integrate(dnorm, -1, 1, mean=0, sd=1)
P <- function(x, mean, sd)
{
  variance <- sd^2
  exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance)
}


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

# > apply(d[,-1],BY.COL<-2,sum)
# AGE_1 AGE_2 AGE_3 AGE_4 AGE_5 TOTAL
#    55   243   156    47    22   523

weight <- c(55, 243, 156, 47, 22)  # Number of fish
weight <- weight / sum(weight)

mean   <- c(23.32727, 33.09053, 41.27244, 51.24468, 61.31818)

sigma  <- c(2.436638, 2.997594, 4.274585, 5.077521, 7.070303)

Delta <- 0.2
x <- seq(0,80, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y4 <- weight[4]*P(x, mean[4], sigma[4])
y5 <- weight[5]*P(x, mean[5], sigma[5])
y <- y1 + y2 + y3 + y4 + y5

par(mfrow=c(1,1))

plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")
lines(x,y4, col="magenta")
lines(x,y5, col="cyan")


derivative1 <- Deriv1(x,y)
derivative2 <- Deriv2(x,y)

par(mfrow=c(3,1))
plot(x,y, type="l", lwd=3,
  main="Heming Lake Pike: Distribution by Age Groups",
  xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")

lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")
lines(x,y4, col="magenta")
lines(x,y5, col="cyan")

##### 1st Derivative

plot(derivative1$x,derivative1$y, type="l",
  main="1st Derivative", xlab="x", ylab="y'")
abline(h=0, lty="dotted")

peaks.Deriv1   <- peaks( derivative1$y, span=3)
valleys.Deriv1 <- peaks(-derivative1$y, span=3)

points( derivative1$x[peaks.Deriv1], derivative1$y[peaks.Deriv1],
        pch=19, col="red")
points( derivative1$x[valleys.Deriv1], derivative1$y[valleys.Deriv1],
        pch=19, col="blue")

# Approximate location of peak and valley
derivative1$x[peaks.Deriv1]
derivative1$x[valleys.Deriv1]

s.approx <- (derivative1$x[valleys.Deriv1] - derivative1$x[peaks.Deriv1])/2
s.approx

# Approximate standard deviation
(derivative1$x[valleys.Deriv1] -  derivative1$x[peaks.Deriv1])/2

##### 2nd Derivative

plot(derivative2$x,derivative2$y, type="l",
  main="2nd Derivative", xlab="x", ylab="y''")
abline(h=0, lty="dotted")

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


# Attempt non-linear curve fit to extract the parameters
set.seed(17)
y <- y + rnorm(length(y), 1E-5, 1E-4)

fit.pike <- nls(y ~
                (a/b)*exp(-(x-c)^2/(2*b^2)) +
                (d/e)*exp(-(x-f)^2/(2*e^2)) +
                (g/h)*exp(-(x-i)^2/(2*h^2)) +
                (j/k)*exp(-(x-l)^2/(2*k^2)),
                start=list(a=(1/sqrt(2*pi)) / s.approx[1], b=s.approx[1], c=mu.approx[1],
                           d=(1/sqrt(2*pi)) / s.approx[2], e=s.approx[2], f=mu.approx[2],
                           g=(1/sqrt(2*pi)) / s.approx[3], h=s.approx[3], i=mu.approx[3],
                           j=(1/sqrt(2*pi)) / s.approx[4], k=s.approx[4], l=mu.approx[4] ),
        control=nls.control(tol=1E-5, minFactor=1/1024),
        trace=TRUE)

# Means (mu values)
coef(fit.pike)[3*1:4]

# Standard deviations (s values)


#dev.off()
