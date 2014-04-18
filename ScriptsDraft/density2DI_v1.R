<<<<<<< HEAD
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


y = c(den1.dt,den2.dt)
plot(density(y))


#Sample data with unbalanced proportion

dat.2 <- data.frame(dens = c(rnorm(1000), rnorm(10, 10, 5))  , lines = c(rep("a",1000), rep( "b", 10)))
densityplot(~dens, data=dat.2, groups = lines,  plot.points = FALSE, ref = TRUE, auto.key = list(space = "right"))

plot(density(dat.2$dens), main = "Density of mixture of two data sets")
plot(density(dat.2$dens[which(dat.2$lines=="a")]), main = "Density of 1000 norm sim data")
plot(density(dat.2$dens[-which(dat.2$lines=="a")]), main = "Density of 10 norm (10,5) sim data")


#Plot.


library(ggplot2)
ggplot(dat.2, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)

library(lattice)
densityplot(~dens,data=dat.2, groups = lines, plot.points = FALSE, ref = TRUE, auto.key = list(space = "right"))

##========================================================
#	On a real D.I. value data
##=========================================================


f_IN <-  "x:/myGit/mixturemodel/data/128110.csv"
dt <- read.csv (f_IN)
title = paste ("sample_", sub(".csv", "", f), sep="")
plot(density(as.vector(dt$DNA_Index)), ylab = "Density", xlab = "DNA Index value", main= title )
=======
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
>>>>>>> 5c56ee008cd8af9fe910e2f6e2bb19f91da2ab8e
