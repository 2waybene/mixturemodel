source("x:/R-project/customPackages/dataManipTools.R")
source("x://R-project/customPackages/arraySeqTools.R")
source("x:/R-project/customPackages/plotTools.R")
source("X:/project2012/ASE/R-scripts/functions-ASE.R")



library(gplots)
library(Rlab)
library(affy)


scriptDir <- "X:/myGit/mixturemodel/workingScripts/"
workingDir  <- "X:/myGit/mixturemodel/workingDir/"
dataDir <- "X:/myGit/mixturemodel/data/dt_12032013/"
dataDir <- "X:/myGit/mixturemodel/data/dt_01232014/"

dataDir <- "/Users/li11/myGit/mixturemodel/data/dt_01232014/"

##=====Start from here..
##  Stat
##=================================================================
dataDir <- "/Users/li11/myGit/mixturemodel/data/dt_01232014/OLK"
workingDir <- "/Users/li11/myGit/mixturemodel/workingDir/"
setwd(workingDir)

fn <- list.files (dataDir, pattern = ".csv")
numOfFile <- length(fn)
i =1 

fn[i]


	dt <- read.csv (paste (dataDir, fn[i], sep="/"))
	#plotDensity(dt$DNA_Index, ylab = "Density", xlab = "DNA Index value", main= title )
plot(density(as.numeric(as.vector(dt$DNA_Index))))









