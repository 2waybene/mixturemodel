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
densityplot(~dens, data=dat.2, groups = lines,  plot.points = FALSE, ref = TRUE, auto.key = list(space = "right"), main = "Two populations in 10:1 ratio")

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


f_IN <-  paste (root, "myGit/mixturemodel/data/128110.csv", sep ="")
dt <- read.csv (f_IN)
title = paste ("sample_", sub(".csv", "", f), sep="")
title = "Sample_128110 density"
plot(density(as.vector(dt$DNA_Index)), ylab = "Density", xlab = "DNA Index value", main= title )


get.den <- density(as.vector(dt$DNA_Index))
plot(get.den, ylab = "Density", xlab = "DNA Index value", main= title )

peak.quick <- function (x, y){
  
  return(x[which(diff(sign(diff(y)))==-2)])
}



str(get.den)
peaks <- peak.quick (get.den$x, get.den$y)
str(peaks)

plot(get.den$x, get.den$y)
get.den$x < peaks[1]

##===================================
# Let's see how far can I do in R
##===================================

dt.raw <- as.vector (dt$DNA_Index)
summary(dt.raw)
dt.raw - peaks[1]
length(dt.raw)

dt.normed = dt.raw - peaks[1] 
which(dt.normed > 0)
which(dt.normed < 0)

dt.first.left <- dt.normed[which(dt.normed < 0)]
dt.first.right <- -dt.first.left
dt.first <- c(dt.first.left, dt.first.right)
plot(density(dt.first))
plot(density(dt.first + peaks[1]))
sim.dt.first <- rnorm(1778, peaks[1], (sd(dt.first+peaks[1])))
plot(density(sim.dt.first), bw=0.01613)

summary(sim.dt.first)
which(dt.raw == (dt.first + peaks[1]))
summary(dt.first + peaks[1])
mean(dt.first+peaks[1])
sd(dt.first + peaks[1])

summary(dt.raw)

peak.quick(density(sim.dt.first)$x, density(sim.dt.first)$y)


##  Fact about first part
first.mean <- mean(dt.first+peaks[1])
first.sd   <- sd(dt.first + peaks[1])
first.den  <- density(dt.first + peaks[1])
first.sim.den <- density(sim.dt.first)



##  standardize den$y, so that the integral will equal 1
sum(first.sim.den$y)
tem <- first.sim.den
tem$y <- first.sim.den$y/365.0299
plot(tem)
str(tem)


((tem$y[256]+tem$y[257])/2)*1778
 