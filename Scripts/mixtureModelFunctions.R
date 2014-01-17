##===============================================================
#	File: mixtureModelFunctions.R
#	Author: Jianying Li
# 	Initial coded by: efg, 23 Aug 2006 at Stower Institutes
##=============================================================== 		

#	P (x, mean, sd)
#		x 	- given any "x"
#		mean 	- normdist mean
#		sd	- normdist standard deviation
#		return the "density probability"

# 	Here P is the same as dnorm (probability density function for normal
#	distribution), but other functions could be tried here.
# 	integrate(P, -1, 1, mean=0, sd=1) is same as integrate(dnorm, -1, 1, mean=0, sd=1)
# 	integrate(P,-5,5, mean=0, sd=1)   # should be close to 1.0


#	peaks(series,span=3)
#		series	- a vector of number
#		span 		- defauty 3
#		return (where a peak is found)

# 	Find "peaks" in array.
# 	R equivalent of Splus peaks() function.
# 	http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33097.html
# 	(see efg's posting to R-Help on 8 Feb 2007 about problem with ties.)
#
#	 peaks(c(1,4,4,1,6,1,5,1,1),3)
# 	[1] FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE

##========Implementation starts here....===============================

P <- function(x, mean, sd)  #-- Test DONE --#
{
  variance <- sd^2
  exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance)
}



peaks <- function(series,span=3)   #-- Test DONE --#
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
Deriv1 <- function(x,y)  #-- Test DONE --#
{
  y.prime <- diff(y) / diff(x)
  x.prime <- x[-length(x)] + diff(x)/2
  list(x = x.prime,
       y = y.prime)
}

# "Centered" 2nd-derivative. Spacing of x-points assumed to be uniform.
Deriv2 <- function(x,y)  #-- Test DONE --#
{
  h <- x[2] - x[1]
  Range <- 2:(length(x)-1)  # Drop first and last points
  list(x = x[Range],
       y = (y[Range+1] - 2*y[Range] + y[Range-1]) / h^2)
}

##==================================##
##	Peak functions here		##
##==================================##


##	peak.quick	##-- Test DONE --#


peak.quick <- function (x, y){

	return(x[which(diff(sign(diff(y)))==-2)])
}


##==============================================
# First derivative from splinefun 
# Same as Deriv1
# The following functions are yet to be tested
##==============================================

firstDerivative = function (x, y ) {

  # fit a smooth spline, and return a function describing it
  fx.spline <- splinefun(x, y)
  # get the negative of 1st derivative
  fx.neg1st <- -fx.spline(x, deriv = 1)
  return(fx.neg1st)
} # end of function



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


