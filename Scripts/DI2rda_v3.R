##==============================================================
##  FileName: DI2rda_v3.R
##  Author: Jianying Li
##  Comment: this is "n" version of the script with goals
##           (1) Determine how many populations to save
##           (2) Thresholding for three populations
##              a. Normal family: mean < 1.2
##              b. Mitotic family: mean < 2.3
##              c. Aneuploidy family: mean > 2.3
##           (3) Upto two round cleaning for normal family
##           (4) Upto three round cleaning for mitotic family
##           (5) Output an .rda object, which contains:
##              -a. sample ID
##              a. count, mean and std for normal family
##              b. count, mean and std for mitotic family
##              c. D.I. values for aneuploidy family
##		So far, works (v1) perfectly fine with OSCC data!!
##		And, works (v2) perfectly fine with OSCC data!!
##		Okay,works (v3) perfectly fine with normal data
##		Done with data cleaning!!
##=============================================================


##  library and customized functions
library(Rlab)

##  OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- linux
#root <- windows
root <- mac.os

source (paste (root, "myGit/mixturemodel/Scripts/cleaningFuncs.R", sep = ""))
##=============================================
# Three criteria and can be modified
# aneuploidy: > 2.3c --> three populations
# mitotic:    > 1.7c --> two populations
# normal??

aneuThresh = 2.3
mitoThresh = 1.7
mitoMin    = 1.5
normThresh = 0.9 ## Not necessary at some point
normMax    = 1.2 ## For additional round of cleaning

##=============================================
##  Read in data
##=============================================

rawFiles <- list.files (paste(root, "myGit/mixturemodel/data/dt_01232014/Normal/",sep=""), pattern = "csv")
rawFiles


for (i in 1:length(rawFiles))
{


	fileName <- paste("myGit/mixturemodel/data/dt_01232014/Normal/", rawFiles[i], sep ="")
	f_IN <-  paste (root, fileName, sep ="")
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

## determine how many families are we dealing with
	numOfFamily <-  1 # minimun one family

	if (length(which(as.vector(dt$DNA_Index) > aneuThresh)) >= 1)
	{
  		numOfFamily = 3
	}else if (length(which(as.vector(dt$DNA_Index) > mitoThresh)) > 1)
	{
  		numOfFamily = 2
	}

numOfFamily

##===================================================
##  removing the normal family
##  upto two round
##===================================================

	dt.raw  <- ""
	firstDT <- ""
	get.gen <- ""
	peaks   <- c()

	dt.raw <- as.vector (dt$DNA_Index)
	tol.num.of.dt <- length(dt.raw)
	get.den <- density(dt.raw)
	peaks <- peak.quick (get.den$x, get.den$y)
peaks

##==========================================================
##  Determine where to start the first population
##  There could be more small peaks less than 1
##  Try to get the first one peaks > 1 but < 1.2 (normMax)
##==========================================================

	index = 1

	length(which(peaks < 1))

	if (!(is.na((peaks[length(which(peaks<1)) + 1]))) & peaks[length(which(peaks<1)) + 1] < normMax ) 
	{ 
  		index = length(which(peaks<1)) + 1
	}else { index = length(which(peaks<1)) }

index

##============================================
##  clean starts here with first population
##============================================

	firstDT <- getPopWIndex (dt.raw, index)

	##  Save first population dt
	FP_dt_primary <- firstDT + peaks[index]
	dt.cleaned <- cleanFirstPop(peaks[index], firstDT, dt.raw)

	if (length(firstDT) >= length(dt.raw) & numOfFamily <=2)
	{
		SP_dt_primary <- dt.raw[c(which(dt.raw > normMax))]
		FP_mean  <- mean(FP_dt_primary)
		FP_std   <- sd(FP_dt_primary)
		FP_count <- length(FP_dt_primary)
		FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
		cleanedSample <- c(cleanedSample, FP)
	}else{

##================================
##  Second round if ever needed
##================================
		if (length (get.den <- tryDensity (dt.cleaned)) >=1 )	
		{
			dt.raw  <- dt.cleaned
			firstDT <- ""
			get.gen <- ""
			peaks   <- c()

			index = 1
			firstDT <- getPopWIndex (dt.raw, index) ##FIXME, yep breaks at sample 88!!!!
			get.den <- density(dt.raw)
			peaks <- peak.quick (get.den$x, get.den$y)
peaks

##=========================================
##  Follow the same protocol, but just 
##  carry out one more cleaning cycle
##  if there is any peaks less than 1.2(normMax)
##===========================================

##Need to add the "cleaned back to population one"!!

			dt.clean.return = list()
			if (peaks[1] < normMax )
			{
  				dt.clean.return <- followUpClean (peaks[1], firstDT, dt.raw)
	    			if (length(dt.clean.return$dtFiltered) > 1) #FIXME, there was a bug
  				{
    					FP_dt_primary <- c(FP_dt_primary, dt.clean.return$dtFiltered)
  				}
  				dt.1pop.cleaned <- dt.clean.return$dtRetain
			}else{
  				dt.1pop.cleaned <- dt.cleaned
			}

##===================================
##  Storing the cleaning results
##===================================

			FP_mean  <- mean(FP_dt_primary)
			FP_std   <- sd(FP_dt_primary)
			FP_count <- length(FP_dt_primary)
			FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
			cleanedSample <- c(cleanedSample, FP)
		} else
		{
			FP_mean  <- mean(FP_dt_primary)
			FP_std   <- sd(FP_dt_primary)
			FP_count <- length(FP_dt_primary)
			FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
			cleanedSample <- c(cleanedSample, FP)
			dt.1pop.cleaned <- dt.cleaned 
		}

	}

##===========================================
##  No need to clean the second population
##===========================================
		SP_mean  <- mean(SP_dt_primary)
		SP_std   <- sd(SP_dt_primary)
		SP_count <- length(SP_dt_primary)
		SP <- list ("SP_mean" = SP_mean, "SP_std" = SP_std, "SP_count" = SP_count)
		cleanedSample <- c(cleanedSample, SP)

		aneup.pop <- ""
		aneu <- list ("AneuLeft" = aneup.pop)
		cleanedSample <- c(cleanedSample, aneu)
		cleanedSample

##==========================
##  Saving the results
##==========================
	storage.dir <- paste (root, "myGit/mixturemodel/cleanedData/Normal/", sep = "")
	file2save <- paste (storage.dir, "cleaned_", cleanedSample$sample, ".rda", sep="")
	save (cleanedSample, file = file2save)

}


