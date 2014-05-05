##	learningDensity.R

library(lattice)
dat <- data.frame(dens = c(rnorm(100), rnorm(100, 10, 5))
                   , lines = rep(c("a", "b"), each = 100))

densityplot(~dens,data=dat,groups = lines,
            plot.points = FALSE, ref = TRUE, 
            auto.key = list(space = "right"))

#	File: simulate_figures




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
            plot.points = FALSE , ref = TRUE ,
		main = "Two populations in 1:1 ratio",
		xlab = "X vaules"
)

leg = c("A population", "B population")
legend ("topright", leg,  col = c("blue", "red"), lty = 1)




#library(ggplot2)
#ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)

##	Assuming mixture case

#plot(density(dat$dens))

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
mtext ("Two population in 100:1 ratio")


dat.2 <- data.frame(dens = c(rnorm(1000, 1,0.5), rnorm(10, 2, 0.5))  , lines = c(rep("a",1000), rep( "b", 10)))
densityplot(~dens, data=dat.2, groups = lines,  plot.points = FALSE, ref = TRUE, auto.key = list(space = "right"))
plot(density(dat.2$dens), main = "Density of mixture of two data sets 100:1 ratio")
mtext ("Equal variance, mean2 is just two times of mean1")


plot(density(dat.2$dens[which(dat.2$lines=="a")]), main = "Density of 1000 norm sim data")
plot(density(dat.2$dens[-which(dat.2$lines=="a")]), main = "Density of 10 norm (10,5) sim data")


#Plot.


library(ggplot2)
ggplot(dat.2, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)

library(lattice)
densityplot(~dens,data=dat.2, groups = lines, plot.points = FALSE, ref = TRUE, auto.key = list(space = "right"))


