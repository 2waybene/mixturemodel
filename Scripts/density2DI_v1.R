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

##===============================
##	Can we find the peak
##===============================
which(dat$lines=="a")
den1.dt <- dat$dens[which(dat$lines=="a")]
den2.dt <- dat$dens[which(dat$lines=="b")]

den1 <- density(den1.dt)
den2 <- density(den2.dt)


y = den1$y + den2$y

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
