##==============================================================
##	File: testing_functions.R
##	Author: Jianying Li
##	Comment: This script is to test customized functions
##		mixtureModelFunctions.R
##		peakFunctions.R
##===============================================================
source("x:/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")

##========================
##	Getting data ready
##========================
Delta <- 0.2
x <- seq(0,80, by=Delta)
mean   <- c(23.32727, 33.09053, 41.27244, 51.24468, 61.31818)
sigma  <- c(2.436638, 2.997594, 4.274585, 5.077521, 7.070303)
weight <- c(55, 243, 156, 47, 22)  # Number of fish
weight <- weight / sum(weight)

##===Testing function "P" here=================
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3 

par(mfrow=c(1,1))

plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")

##==========================================
##	Testing function "Deriv1 & Deriv2"
##	First/second derivative
##==========================================
derivative1 <- Deriv1(x,y)
derivative2 <- Deriv2(x,y)

#	DV1 <- firstDerivative(x,y)


par(mfrow=c(3,1))
plot(x,y, type="l", lwd=3,
  main="Heming Lake Pike: First three groups only",
  xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")

lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")

##### 1st Derivative

plot(derivative1$x,derivative1$y, type="l",
  main="1st Derivative", xlab="x", ylab="y'")
abline(h=0, lty="dotted")

##===Testing function "peaks" here=================
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


##==============================================
##	Testing peak calling functions
##==============================================
##	Three (really two) peaks data
mean <- c(1.001, 2.002, 2.300)
sigma <- c(0.19, 0.25,0.5) 
weight <- c(0.893, 0.092, 0.05)


Delta <- 0.01
x <- seq(0,7, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3

plot(x,y) 	## Only TWO visible peaks

##======================##
##	A better plot	##
##======================##

par(mfrow=c(1,1))
plot(x,y, type="l", lwd=3,
     main="Simulated D.I. values: three-category groups",
     xlab="D.I. value", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")

lgd = c("Mixture", "Normal, mean = 1.001", "Mitotic,  mean = 2.002", "Aneuploid, mean = 2.300")
legend ("topright", lgd, text.col = c("black", "red", "green", "blue"))


##	Calling peaks...

peakx <- peak.quick(x,y)
peakx

peakx <- peak.quick(x,y1)
peakx #found the obvious peaks, but miss the trivil one(s)!!



