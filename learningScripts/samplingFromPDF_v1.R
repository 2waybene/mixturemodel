##=========================================================================================
#	samplingFromPDF_v1.R
#	very interesting sampling from pdf
#	http://blog.quantitations.com/tutorial/2012/11/20/sampling-from-an-arbitrary-density/
##=========================================================================================
invcdf <- function(u, m) {
    return(sqrt(m^2/(1 - (1 - m^2) * u)))
}

sample1 <- sapply(runif(100), invcdf, m = .5)
plot(density(sample1), main = "Sample Density using invcdf Function")

sample1


endsign <- function(f, sign = 1) {
    b <- sign
    while (sign * f(b) < 0) b <- 10 * b
    return(b)
}

samplepdf <- function(n, pdf, ..., spdf.lower = -Inf, spdf.upper = Inf) {
    vpdf <- function(v) sapply(v, pdf, ...)  # vectorize
    cdf <- function(x) integrate(vpdf, spdf.lower, x)$value
    invcdf <- function(u) {
        subcdf <- function(t) cdf(t) - u
        if (spdf.lower == -Inf) 
            spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf) 
            spdf.upper <- endsign(subcdf)
        return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
    }
    sapply(runif(n), invcdf)
}


h <- function(t, m) {
    if (t >= m & t <= 1) 
        return(2 * m^2/(1 - m^2)/t^3)
    return(0)
}

sample2 <- samplepdf(100, h, m = .5)
plot(density(sample2), main = "Sample Density using samplepdf Function")




ntime <- function(n, f, ...) {
    return(system.time(f(n, ...))[3])
}

n <- 100 * 1:10
times <- sapply(n, ntime, f = samplepdf, pdf = h)
fit <- lm(times ~ n)$coeff
fit <- signif(fit, 2)
print(paste("The time needed (seconds) to sample from this density is about ", 
    fit[1], " + ", fit[2], "*n.", sep=""))


## [1] "The time needed (seconds) to sample from this density is about 3.2 + 0.14*n."


plot(n, times, ylab = "seconds", xlab = "n", main = "Time Required", col = 3)
abline(fit, col = 2)


