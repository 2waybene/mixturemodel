##	need to define root: mac, linux or windows

dir <- "/myGit/mixturemodel/oldScripts/"
dir <- paste (root, dir, sep="")
source (paste(dir, "HemingLakePike.R", sep=""))
y <- y1 + y2 + y3 

par(mfrow=c(1,1))

plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")


derivative1 <- Deriv1(x,y)
derivative2 <- Deriv2(x,y)

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

