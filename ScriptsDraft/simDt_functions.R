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


