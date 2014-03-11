<<<<<<< HEAD
##=============================================================
##	Script: seqDataFiler.R
##	Author: Jianying Li
##	Comment: originally coded as density2DI_vX.R 
##		   the script should:
##		(1) sequentially stripping off "noise"
##		(2) storing the summary statistics
##		(3) mixing the distributions
##=============================================================

##	OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

##	package needed
root <- windows



##========================================================
#	Read in a real D.I. value data
##=========================================================


f_IN <-  "x:/myGit/mixturemodel/data/dt_01232014/9559229.csv"
dt <- read.csv (f_IN)
plot(density(as.vector(dt$DNA_Index)), ylab = "Density", xlab = "DNA Index value" )
str(dt$DNA_Index)



##=========================================================
##	Getting the peaks
##==========================================================
get.den <- density(as.vector(dt$DNA_Index))
plot(get.den, ylab = "Density", xlab = "DNA Index value")

peak.quick <- function (x, y){
  
  return(x[which(diff(sign(diff(y)))==-2)])
}


str(get.den)
peaks <- peak.quick (get.den$x, get.den$y)
str(peaks)


##=========================================================
##	Stripping the normal
##==========================================================

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



summary(dt.first + peaks[1])
mean(dt.first+peaks[1])
sd(dt.first + peaks[1])



##  Fact about first part
first.mean <- mean(dt.first+peaks[1])
first.sd   <- sd(dt.first + peaks[1])
first.den  <- density(dt.first + peaks[1])


tem <- first.den
tem$y <- first.den$y/sum(first.den$y)
plot(tem)
str(tem)


##  standardize den$y, so that the integral will equal 1


##  Filter starts here...

##	Retain data on the right of the first peak
##			SAME AS
##	Remove data on the left  of the fist peak

dt.raw.flt.1 <- dt.raw[which(dt.raw >=  peaks[1])]
str(dt.raw.flt.1)

##	Retain data on the right of max of the first population
max(dt.first+peaks[1])
dt.raw.02 <- dt.raw.flt.1[which(dt.raw.flt.1 >=max(dt.first+peaks[1]))]
str(dt.raw.02)


##	Data fall between the right of the first peak and the left of max of the first population
dt.raw.flt.2 <- dt.raw.flt.1[-which(dt.raw.flt.1 >=max(dt.first+peaks[1]))]
summary(dt.raw.flt.2)
str(dt.raw.flt.2)




##	Now, remove the data according to the "estimated proportion"
##	Between two adjacent populations



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


num2salvage <- sample (c(1:length(dt2filter)), ceiling(adjust)) 	#FIXME: Manully fixting
dt2filter <- dt2filter[-num2salvage]

##	Need to get this manual adjusting number ??
num2salvage <- sample (c(1:length(dt2filter)), 138)			#FIXME: Manully fixting!!!
dt2filter <- dt2filter[-num2salvage]
str(dt2filter)

dt.raw.02 <- c(dt.raw.02, dt2filter)
str(dt.raw.02)
plot(density(dt.raw.02), main ="Density after removing the first population")

dt.second.pop <- dt.raw.02
str(dt.second.pop)





##=========================================================
##	Stripping the mitotics
##==========================================================



