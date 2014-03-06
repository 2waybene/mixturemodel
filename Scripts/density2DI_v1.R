##=============================================================
##	Script: density2DI_v1.R
##	Author: Jianying Li
##	Comment: finalized from density, mixture model, and 
##		   peak finding and will be applied to D.I.
##	
##=============================================================

##	OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

##	package needed
root <- windows



##	simulated data and density plot

library(lattice)
dat <- data.frame(dens = c(rnorm(100), rnorm(100, 10, 5))
                   , lines = rep(c("a", "b"), each = 100))
densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE, 
            auto.key = list(space = "right"))

library(ggplot2)
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)


##	Assuming mixture case

plot(density(dat$dens))

mix.den <- density(dat$dens)
str(mix.den)

##===============================
##	Can we find the peak
##=========================================
library(Rlab) #dataset precip is masked!!!

which(dat$lines=="a")
den1.dt <- dat$dens[which(dat$lines=="a")]
den2.dt <- dat$dens[which(dat$lines=="b")]

summary(den1.dt)
stats(den1.dt)
summary(den2.dt)
stats(den2.dt)

mix.again <- c(den1.dt, den2.dt)
str(mix.again)

plot(density(mix.again))

den <- density(mix.again)
str(den)

peak.quick(den$x, den$y)
den1 <- density(den1.dt)
stats(den1)
den2 <- density(den2.dt)
stats(den2)

y = den1$y + den2$y

plot(density(y))
##	D.I. values



#Sample data
dat <- data.frame(dens = c(rnorm(100), rnorm(100, 10, 5))
                   , lines = rep(c("a", "b"), each = 100))
#Plot.

windows()
plot(density(dat$dens, group = dat$lines))

str(dat)

dat$dens
summary(dat$lines)



library(affy)
plotDensity(dat$dens)
