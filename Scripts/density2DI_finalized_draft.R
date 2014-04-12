##=============================================================
##	Script: density2DI_v4.R
##	Author: Jianying Li
##	Comment: finalized from density, mixture model, and 
##		   peak finding and will be applied to D.I.
##	
##=============================================================
library(Rlab)

##	OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

##	package needed



##=============================================
##  Read in data
##=============================================
root <- windows
root <- mac.os
fileName <- "myGit/mixturemodel/data/128110.csv"

##========================================================
#  On a real D.I. value data
##=========================================================
source (paste (root, "myGit/mixturemodel/Scripts/cleanFuncs.R", sep = ""))


f_IN <-  paste (root, fileName, sep ="")


#f_IN <-  paste (root, "myGit/mixturemodel/data/128110.csv", sep ="")
#f_IN <-  paste (root, "myGit/mixturemodel/data/dt_01232014/9559572.csv", sep ="")

nameSplit <- strsplit(f_IN, "/")[[1]]
sampleName <- nameSplit[length(nameSplit)]
sampleName <- sub(".csv", "", sampleName)
sampleName

cleanedSample <-  list("sample" = sampleName)
cleanedSample


##============================
# read in the raw D.I. value
##============================

dt <- read.csv (f_IN)
get.den <- density(as.vector(dt$DNA_Index))

peak.quick <- function (x, y){
  
  return(x[which(diff(sign(diff(y)))==-2)])
}


peaks <- peak.quick (get.den$x, get.den$y)
str(peaks)


##===================================
# getting the first population stats
##===================================

dt.raw <- as.vector (dt$DNA_Index)
summary(dt.raw)
stats(dt.raw)
#dt.raw - peaks[1]
totalCount <- length(dt.raw)

dt.normed = dt.raw - peaks[1] 
dt.first.left <- dt.normed[which(dt.normed < 0)]
dt.first.right <- -dt.first.left
dt.first <- c(dt.first.left, dt.first.right)

#plot(density(dt.first))
#plot(density(dt.first + peaks[1]), main = "Density from the left most population")


FP_mean <- mean(dt.first+peaks[1])
FP_std <- sd(dt.first + peaks[1])
FP_count <- length(dt.first + peaks[1])

FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)

cleanedSample <- c(cleanedSample, FP)
cleanedSample

##======================================================
##  Now, let's work on removing the first population
##======================================================

##	clean with function


dt.cleaned <- cleanFirstPop ( peaks[1],   dt.first,   dt.raw  ) 
plot(density(dt.cleaned), main ="Density after removing the first population")


firstDT <- getFirstPop(dt.cleaned)
plot(density(firstDT))
peaks <- peak.quick(density(dt.cleaned)$x, density(dt.cleaned)$y)

dt.another.clean <- cleanFirstPop(peaks[1], firstDT, dt.cleaned)
plot(density(dt.another.clean))


##=============================================
# Take care of the residual 
##=============================================



##=============================================
# result from first filtering
# dt.second.pop <- dt.raw.02
##=============================================

get.den <- density(dt.second.pop)
plot(get.den, ylab = "Density", xlab = "DNA Index value", main= title )


get.den <- density(dt.raw.02)

str(get.den)
peaks <- peak.quick (get.den$x, get.den$y)
str(peaks)


##===================================
# Let's see how far can I do in R
##===================================

dt.raw <- dt.second.pop
summary(dt.raw)
dt.raw - peaks[2]
length(dt.raw)

dt.normed = dt.raw - peaks[2] 
which(dt.normed > 0)
which(dt.normed < 0)

dt.first.left <- dt.normed[which(dt.normed < 0)]
dt.first.right <- -dt.first.left
dt.first <- c(dt.first.left, dt.first.right)
plot(density(dt.first))
plot(density(dt.first + peaks[2]))
str(dt.first)
sim.dt.first <- rnorm(1308, peaks[2], (sd(dt.first+peaks[2])))
plot(density(sim.dt.first), bw=0.01613)

summary(sim.dt.first)
which(dt.raw == (dt.first + peaks[2]))
summary(dt.first + peaks[2])
mean(dt.first+peaks[2])
sd(dt.first + peaks[2])

summary(dt.raw)

peak.quick(density(sim.dt.first)$x, density(sim.dt.first)$y)


##  Fact about first part
first.mean <- mean(dt.first+peaks[2])
first.sd   <- sd(dt.first + peaks[2])
first.den  <- density(dt.first + peaks[2])
first.sim.den <- density(sim.dt.first)
save (first.mean, first.sd , file = "sample_128110_stat_pop2.rda")

##======================================================
##  Now, let's work on removing the first population
##======================================================

##  standardize den$y, so that the integral will equal 1
sum(first.sim.den$y)
first.den <- density(dt.first + peaks[2])
str(first.den)
sum(first.den$y)


tem <- first.den
tem$y <- first.den$y/sum(first.den$y)
plot(tem)
str(tem)


##  Filter starts here...

##	Retain data on the right of the first peak
##			SAME AS
##	Remove data on the left  of the fist peak

dt.raw.flt.1 <- dt.raw[which(dt.raw >=  peaks[2])]
str(dt.raw.flt.1)


##	Retain data on the right of max of the first population
max(dt.first+peaks[2])
dt.raw.02 <- dt.raw.flt.1[which(dt.raw.flt.1 >=max(dt.first+peaks[2]))]
str(dt.raw.02)


##	Data fall between the right of the first peak and the left of max of the first population
dt.raw.flt.2 <- dt.raw.flt.1[-which(dt.raw.flt.1 >=max(dt.first+peaks[2]))]
summary(dt.raw.flt.2)
str(dt.raw.flt.2)

##	Now, remove the data according to the "estimated proportion"
##	Between two adjacent populations

plot(density(tem$x))
str(tem$x)


adjust = 0; 
dt2filter <- dt.raw.flt.2
str(dt2filter)
for (i in 1:256)
{
temp = 0
	l.bound <- i + 255
	h.bound <- i + 256
	num.of.data <- ((tem$y[l.bound] + tem$y[h.bound])/2)*(length(dt.first))
	candidate <- which(dt2filter > tem$x[l.bound] & dt2filter  < tem$x[h.bound])

	if (length(candidate) >=1)
  	{
    		if (length(candidate) > floor(num.of.data))
    		{
			temp = num.of.data - floor(num.of.data)
      		data2exclude <- sample(candidate, floor(num.of.data))
      		if (length(data2exclude) >=1 )
      		{
        			dt2filter <- dt2filter[-data2exclude]	
				adjust = adjust + temp
      		}
		}else{
      		dt2filter <- dt2filter[-candidate]
    		}
  	}
}

str(dt2filter)
adjust

num2salvage <- sample (c(1:length(dt2filter)), ceiling(adjust)) 	#FIXME: Manully fixting
#num2salvage <- sample (c(1:length(dt2filter)), 138)			#FIXME: Manully fixting!!!
dt2filter <- dt2filter[-num2salvage]

num2salvage <- sample (c(1:length(dt2filter)), length(which(dt2filter < (first.mean + first.sd))))
dt2filter <- dt2filter[-num2salvage]
str(dt2filter)
plot(density(dt2filter))

dt.raw.02 <- c(dt.raw.02, dt2filter)
str(dt.raw.02)
summary(dt.raw.02)
plot(density(dt.raw.02), main = "After removing both normal and replicating populations")
stats(dt.raw.02)

get.den <- density(dt.raw.02)
peaks <- peak.quick (get.den$x, get.den$y)
str(peaks)


dt2work <- dt.raw.02

