##=============================================================
##	Script: 	getDensity4all.R
##	Author: 	Jianying Li
##	Comment:	Will loop through all samples and produce  
##		   	density plot for each sample
##	
##=============================================================
library(affy)

##	OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

##	package needed
root <- windows



dataDir <- paste (root, "myGit/mixturemodel/data/dt_01232014/", sep="")
workingDir <- paste (root, "myGit/mixturemodel/workingDir/", sep="")
setwd(workingDir)

fn <- list.files (dataDir, pattern = ".csv")
numOfFile <- length(fn)

for (i in 1:numOfFile)
{
	title = paste ("sample_", sub(".csv", "", fn[i]), sep="")
	outFile = paste (title, ".jpeg", sep="")
	jpeg (outFile)
	dt <- read.csv (paste (dataDir, fn[i], sep=""))
	plotDensity(as.data.frame(dt$DNA_Index), ylab = "Density", xlab = "DNA Index value", main= title )
	dev.off()
}









