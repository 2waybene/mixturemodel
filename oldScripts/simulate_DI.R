##====================================================
##	File:   simulate_DI.R
##	Author: Jianying Li
##	Comment: Learning mixture model for Aneuploidy
##		   study for OSCC
##====================================================


##	Technical hurdles
##	Unbalanced data distribution ??



source("x:/myGit/mixturemodel/ScriptsDraft/mixtureModelFunctions.R")
#source("/Users/li11/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")



##=====Sample data
macFileDir <- "/Users/li11/myGit/mixturemodel/data/"
winFileDir <- "x:/myGit/mixturemodel/data/"

f <- "fmd_DT_raw.txt"
f <- "oscc-olk1_parsed.txt"
f_IN <- paste (file.dir, f , sep="")
f_IN <- paste (macFileDir, f , sep="")
f_IN <- paste (winFileDir, f , sep="")
##=====End of settings

dt <- read.table(f_IN, header= F, sep = "\t")

str(dt)
temp.cv <- sd(dt$V1)/mean(dt$V1)
temp.cv

##=====End

##======================================================================
##  Based on our prior knowledge, we can focus on 
##  three clusters of cell populations:
##  Cluster one, with mean DI = 1.001
##  Cluster two, with mean DI = 2.001
##  Cluster three, with minimum DI > 2.300, assigning 2.300 as mean
##  Now, let's assume that C.V. is same 
##======================================================================

mean <- c(1.001, 2.002, 2.300)
mean <- c(1.001, 2.002, 3.600)
mean <- c(0.936, 1.692, 3.600)
sigma <- c(0.19, 0.25,0.5) 
sigma <- mean*temp.cv

weight <- c(0.893, 0.092, 0.05)
weight <- c(0.425, 0.425, 0.05)


Delta <- 0.01
x <- seq(0,7, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3


#integrate(P, 0,7, mean[1], sigma[1])
#integrate(P, 0,7, mean[2], sigma[2])
#integrate(P, 0,7, mean[3], sigma[3])

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

derivative1 <- Deriv1(x,y)
derivative2 <- Deriv2(x,y)

par(mfrow=c(3,1))
plot(x,y, type="l", lwd=3,
  main="Simulated D.I. values: three-category groups",
  xlab="D.I. value", ylab="Probability Density")
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


# Peaks
derivative2$x[peaks.Deriv2]
derivative2$y[peaks.Deriv2]


# Valleys
derivative2$x[valleys.Deriv2]
derivative2$y[valleys.Deriv2]



# Attempt non-linear curve fit to extract the parameters
set.seed(17)
y <- y + rnorm(length(y), 1E-5, 1E-4)

fit.pike <- nls(y ~
                (a/b)*exp(-(x-c)^2/(2*b^2)) +
                (d/e)*exp(-(x-f)^2/(2*e^2)),
 
                start=list(a=(1/sqrt(2*pi)) / mean[1], b=mean[1], c=sigma[1],
                           d=(1/sqrt(2*pi)) / mean[2], e=mean[2], f=sigma[2]),
                control=nls.control(tol=1E-5, minFactor=1/1024),
        trace=TRUE)

# Means (mu values)
coef(fit.pike)[3*1:4]

# Standard deviations (s values)


 
##========================================================
## Summary data from 29 samples
##========================================================


f <- "summary_dt.txt";
f_IN <- paste (file.dir, f , sep="")
f_IN <- paste (macFileDir, f , sep="")
f_IN <- paste (winFileDir, f , sep="")
##=====End of settings

dt <- read.table(f_IN, header= T, sep = "\t")

str(dt)

summary(dt$Mean3)



