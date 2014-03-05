##=============================================================
##	File name: extracPeaksOralDt.R
##	Author   : Jianying Li
##	History  : Initial coded, 12/09/2013
##	Comment  : It is important to compare different
##		     peak functions
##=============================================================
source("x:/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")
source("x:/myGit/mixturemodel/Scripts/R-peak-functions.R")
library(affy)
library(Rlab)
winFileDir <- "x:/myGit/mixturemodel/data/"
f <- "128110.csv"
f_IN <- paste (winFileDir, f , sep="")

dt <- read.csv (f_IN)
title = paste ("sample_", sub(".csv", "", f), sep="")
plotDensity(as.data.frame(dt$DNA_Index), ylab = "Density", xlab = "DNA Index value", main= title )

plot(as.data.frame(dt$DNA_Index), ylab = "Density", xlab = "DNA Index value", main= title )



stats(as.data.frame(dt$DNA_Index))
x <- seq(min(as.data.frame(dt$DNA_Index)),max(as.data.frame(dt$DNA_Index)), length.out =  dim(as.data.frame(dt$DNA_Index))[1])
y <- as.data.frame(dt$DNA_Index)

plotDensity(y)
hist(y[[1]])
stats(y[[1]])


##	okay, y[[1]] has 2739
##	Now, let's find some peaks..

p=findPeaks(y[[1]])
p #NO peaks found

p=findPeaks(y[[1]], 1)
peakx <- peak.dt.only(y[[1]],0.1)
which(peakx == TRUE)#NO peaks found


p = peakfinder(y[[1]])
p  #Found peaks but not correct




##==================================##
##	Need peak functions here	#
###==================================##



data <-  y[[1]]

peaks <- peakfinder(data)
hist(data)
plotDensity(y)
sapply(peaks,function(x) abline(v=x,col="red"))
peaks

