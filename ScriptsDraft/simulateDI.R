#	simulating DNA D.I. values
# Need to source ~/myGit/mixturemodel/Scripts/simDt_functions.R
#	normal population
library(Rlab)

##  OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"


##=============================================
##  Read in data
##=============================================
root <- windows
#root <- mac.os

source (paste (root, "/myGit/mixturemodel/Scripts/simDt_functions.R", sep=""))


mean.norm <- c()
for (i in 1:20)
{
	mean.norm[i] <- runif (1, 0.9,1.2)
}


mean.normalP <- mean.norm[sample(1:20,1)]

#	mitotic population
mean.mitotic <- c()
for (i in 1:20)
{
	mean.mitotic[i] <- runif (1, 1.7,2.2)
}
mean.mitoticP <- mean.mitotic[sample(1:20,1)]


#	aneuploidy  population

mean.aneu <- 3.3
st.aneu <-3.5
n = 40

sample.aneu <- rnorm(n, mean = mean.aneu, sd = st.aneu)
sample.aneuP <- sample.aneu[sample.aneu >=0]
plot(density(sample.aneuP))
mean.aneuP <- mean(sample.aneuP)
sigmaA <- sd(sample.aneuP)


##=======================================
##	sample D.I. range
##=======================================
Delta <- 0.001
x <- seq(0,11, by=(11/512))
x <- x[-1]
length(x)

y1 <- P(x, mean.normalP, 0.3)
y2 <- P(x, mean.mitoticP, 0.3)



plot(x, y1)
plot(x, y2)

y3 <- P(x, mean(sample.aneuP), sd(sample.aneuP))
plot(x, y3)



##=============================
##	Different weight
##==============================

##==============================
normWeight <- c(0.66,0.32,0.04)
mitoticWeight <- c(0.48, 0.48,0.04)
OSCCWeight <- c(0.3, 0.2,0.5)
sigmaP <- 0.3
sigmaM <- 0.3
sigmaA <- 3.5

weight <- normWeight
weight <- mitoticWeight
weight <- OSCCWeight


y1 <- weight[1]*P(x, mean.normalP, sigmaP)
y2 <- weight[2]*P(x, mean.mitoticP, sigmaM)
y3 <- weight[3]*P(x, mean.aneuP, sigmaA)

y = y1 + y2
length(y)
length(x)
y = y1 + y2 + y3

plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")



##=============================================
# simulating from density
##=============================================
library(ks)
set.seed(1)


simDt <- rkde(fhat=kde(x=y, h=hpi(y)), n=100, positive=TRUE)
histinfo <- hist(simDt, col="green", freq=F, breaks=50)
histinfo

simDt

##==============================================

abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")



x <- 1:12
# a random permutation
sample(x)

##===============================================
##  A real case
##===============================================
sampleName <- "128123"
storage.dir <- paste (root, "myGit/mixturemodel/cleanedData/OSCC/", sep = "")
file2load <- paste (storage.dir, "cleaned_", sampleName, ".rda", sep="")
load(file2load)


x <- seq(0,9, by=(9/512))
x <- x[-1]
length(x)


##==============================
w.norm <- cleanedSample$FP_count/(cleanedSample$FP_count + cleanedSample$SP_count)*0.9
w.mito <- cleanedSample$SP_count/(cleanedSample$FP_count + cleanedSample$SP_count)*0.9



y1 <- w.norm*P(x, cleanedSample$FP_mean, cleanedSample$FP_std)
y2 <- w.mito*P(x, cleanedSample$SP_mean, cleanedSample$SP_std)

y = y1 + y2
length(y)
length(x)


plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")

x
y

num2sample <- 1000 - length(cleanedSample$AneuLeft)

simDt <- rkde(fhat=kde(x=y, h=hpi(y)), n=num2sample, positive=TRUE)

stats(simDt)


histinfo <- hist(simDt, col="green", freq=F, breaks=15)
histinfo