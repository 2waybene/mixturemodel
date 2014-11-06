##=============================================================
##	Script: 	figure4case_128141.R
##	Author: 	Jianying Li
##	Comment:	Will loop through all samples and produce  
##		   	density plot for each sample
##	
##=============================================================
library(affy)
root <- "X:/"  #this is windows


dataDir <- paste (root, "myGit/mixturemodel/data/dt_01232014/OLK/", sep="")
workingDir <- paste (root, "myGit/mixturemodel/DIdensities/OLK/", sep="")

setwd(workingDir)
getwd()

fn <- list.files (dataDir, pattern = ".csv")
numOfFile <- length(fn)
which(fn == "128141.csv")

dt <- read.csv (paste (dataDir, fn[which(fn == "128141.csv")], sep=""))

plotDensity(as.data.frame(dt$DNA_Index), ylab = "Density", xlab = "DNA Index value", main= title )

##	or

dt$DNA_Index[which(dt$DNA_Index > 2.3)]
get.den <- density(as.vector(dt$DNA_Index))
plot(get.den, xlab = "DNA Index value", main="")

get.den <- density(as.vector(dt.cleaned))
plot(get.den, xlab = "DNA Index value", main="")

##============================
##	EdTAR
##===========================


##======================
##	Second round
##======================

	dt.raw  <- ""
	firstDT <- ""
	get.gen <- ""
	peaks   <- c()

	dt.raw <- dt.cleaned 
	tol.num.of.dt <- length(dt.raw)
	get.den <- density(dt.raw)
	peaks <- peak.quick (get.den$x, get.den$y)
peaks


	numOfFamily <-  1 # minimun one family

	if (length(which(as.vector(dt.raw ) > aneuThresh)) >= 1)
	{
  		numOfFamily = 3
	}else if (length(which(as.vector(dt.raw ) > mitoThresh)) > 1)
	{
  		numOfFamily = 2
	}
numOfFamily



