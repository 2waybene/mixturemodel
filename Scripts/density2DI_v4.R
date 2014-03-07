##=============================================================
##	Script: density2DI_v4.R
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

##========================================================
#	On a real D.I. value data
##=========================================================


f_IN <-  paste (root, "myGit/mixturemodel/data/128110.csv", sep ="")
dt <- read.csv (f_IN)
title = paste ("sample_", sub(".csv", "", f), sep="")
title = "Sample_128110 density"
plot(density(as.vector(dt$DNA_Index)), ylab = "Density", xlab = "DNA Index value", main= title )

# result from first filtering
# dt.second.pop <- dt.raw.02

get.den <- density(dt.second.pop)
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
first.mean <- mean(dt.first+peaks[1])
first.sd   <- sd(dt.first + peaks[1])
first.den  <- density(dt.first + peaks[1])
first.sim.den <- density(sim.dt.first)


##======================================================
##  Now, let's work on removing the first population
##======================================================

##  standardize den$y, so that the integral will equal 1
sum(first.sim.den$y)
first.den <- density(dt.first + peaks[2])
str(first.den)
sum(first.den$y)

#tem <- first.sim.den
#tem$y <- first.sim.den$y/365.0299

tem <- first.den
tem$y <- first.den$y/sum(first.den$y)
plot(tem)
str(tem)


##  Filter starts here...
dt.raw.flt.1 <- dt.raw[which(dt.raw >=  peaks[2])]
str(dt.raw.flt.1)
max(dt.first+peaks[2])
dt.raw.02 <- dt.raw.flt.1[which(dt.raw.flt.1 >=max(dt.first+peaks[2]))]
str(dt.raw.02)

dt.raw.flt.2 <- dt.raw.flt.1[-which(dt.raw.flt.1 >=max(dt.first+peaks[2]))]
summary(dt.raw.flt.2)

plot(density(tem$x))

str(tem$x)

dt2filter <- dt.raw.flt.2
str(dt2filter)
for (i in 1:256)
{
  l.bound <- i + 255
  h.bound <- i + 256
  num.of.data <- (tem$y[l.bound] + tem$y[h.bound])/2*1308
  candidate <- which(dt2filter > tem$x[l.bound] & dt2filter  < tem$x[h.bound])
  if (length(candidate) >=1)
  {
    if (length(candidate) > floor(num.of.data))
    {
      data2exclude <- sample(candidate, floor(num.of.data))
      if (length(data2exclude) >=1 )
      {
        dt2filter <- dt2filter[-data2exclude]
      }
    }else{
      dt2filter <- dt2filter[-candidate]
    }
  }
}

str(dt2filter)

dt.raw.02 <- c(dt.raw.02, dt2filter)
str(dt.raw.02)
plot(density(dt.raw.02))

#dt.second.pop <- dt.raw.02
