##==============================================================
##  FileName: DI2rda_v1 .R
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
normThresh = 1
##=============================================
##  Read in data
##=============================================

rawFiles <- list.files (paste(root, "myGit/mixturemodel/data/dt_01232014/OSCC/",sep=""), pattern = "csv")
rawFiles


for (i in 1:length(rawFiles))
{
  i=6
#print (rawFiles[i])

fileName <- paste("myGit/mixturemodel/data/dt_01232014/OSCC/", rawFiles[i], sep ="")
f_IN <-  paste (root, fileName, sep ="")

nameSplit <- strsplit(f_IN, "/")[[1]]
sampleName <- nameSplit[length(nameSplit)]
sampleName <- sub(".csv", "", sampleName)
sampleName

cleanedSample <-  list("sample" = sampleName)
#cleanedSample


##============================
# read in the raw D.I. value
##============================

dt <- read.csv (f_IN)

## determine how many families are we dealing with
numOfFamily <-  1 # minimun one family

if (length(which(as.vector(dt$DNA_Index) > aneuThresh)) > 1)
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

##===================================================
##  Determine where to start the first population
##  There could be more small peaks less than 1
##  Try to get the first one peaks > 1 but < 1.2
##====================================================

index = 1

length(which(peaks < 1))
if (peaks[length(which(peaks<1)) + 1] < 1.2) 
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

#temp <- followUpClean (peaks[index], firstDT, dt.raw)

#plot(density(dt.cleaned))
#str(dt.cleaned)

##================================
##  Second round if ever needed
##================================
dt.raw  <- dt.cleaned
firstDT <- ""
get.gen <- ""
peaks   <- c()

index = 1


firstDT <- getPopWIndex (dt.raw, index)
#plot(density(firstDT))

get.den <- density(dt.raw)
peaks <- peak.quick (get.den$x, get.den$y)
peaks

##=========================================
##  Follow the same protocol, but just 
##  carry out one more cleaning cycle
##  if there is any peaks less than 1.2
##===========================================

##Need to add the "cleaned back to population one"!!
dt.clean.return = list()

if (peaks[1] < 1.2)
{
  #dt.another.clean <- cleanFirstPop(peaks[1], firstDT, dt.cleaned)
  dt.clean.return <- followUpClean (peaks[1], firstDT, dt.raw)
#  plot(density(dt.another.clean))
#  dt.1pop.cleaned <- dt.another.clean
  #if (length(dt.clean.return$dtFiltered) > 0) #
    if (length(dt.clean.return$dtFiltered) > 1) #FIXME, there was a bug
  {
    FP_dt_primary <- c(FP_dt_primary, dt.clean.return$dtFiltered)
  }
  dt.1pop.cleaned <- dt.clean.return$dtRetain
}else{
  dt.1pop.cleaned <- dt.cleaned
}

FP_mean <- mean(FP_dt_primary)
FP_std <- sd(FP_dt_primary)
FP_count <- tol.num.of.dt  - length(dt.1pop.cleaned)

FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
cleanedSample <- c(cleanedSample, FP)
#cleanedSample


##===========================================
##  Here comes the cleaning for the 
##  second population and store the stats
##===========================================

num.of.DI.left <- length(dt.1pop.cleaned)

get.den <- density(dt.1pop.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
peaks

##  Determine where to start the first population

index = 1

#which(peaks > 1.5)[1]

if (length(which(peaks < 1.5)) >=1 )
{
  index <-which(peaks > 1.5)[1]
 

  #secondDT <- getSecondPop(dt.1pop.cleaned)
  secondDT <- getPopWIndex (dt.1pop.cleaned, index)
#plot(density(secondDT))
##  Save first population stats
SP_dt_primary <- (secondDT + peaks[index])

#  SP_mean <- mean(secondDT + peaks[index])
#  SP_std <- sd(secondDT + peaks[index])

#plot(density(secondDT + peaks[index]))


secondDT.cleaned <- cleanFirstPop(peaks[index],  secondDT, dt.1pop.cleaned)
#str(secondDT.cleaned)
#plot(density(secondDT.cleaned))
}

##==================================
##  Need another round of cleaning
##==================================

get.den <- density(secondDT.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
peaks

##  Determine where to start the first population

index = 0
third_round = 0

if (length(peaks) > 1 &  length(which(peaks < 2)) >= 1)
{
  index =  which(peaks < 2) [length(which(peaks < 2))]
}

#stats (secondDT.cleaned)

#  secondDT.1 <- getFirstPop(secondDT.cleaned)

#secondDT <- getSecondPop(dt.1pop.cleaned)

if (index >=1)
{
  secondDT.1 <- getPopWIndex (secondDT.cleaned, index)

  #plot(density(secondDT.1))


#plot(density(secondDT.1 + peaks[index]))


secondDT.2.cleaned <- cleanFirstPop(peaks[index],  secondDT.1, secondDT.cleaned)
third_round = 1
#str(secondDT.2.cleaned)
#plot(density(secondDT.2.cleaned))
#stats (secondDT.2.cleaned)

}else{
  secondDT.2.cleaned <- secondDT.cleaned
  
}


##  third round??

if (third_round)
{
get.den <- density(secondDT.2.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
#peaks

index = 0

if (length(peaks) > 1 &  length(which(peaks < 1.8)) >= 1)
{
  index =  which(peaks < 1.8) [length(which(peaks < 1.8))]
}

if (index > 0)
{
  secondDT.2.sub <- getFirstPop(secondDT.2.cleaned)
  secondDT.3.cleaned <- cleanFirstPop(peaks[1], secondDT.2.sub , secondDT.2.cleaned)
}else{
  secondDT.3.cleaned <- secondDT.2.cleaned
}
 # plot(density(secondDT.3.cleaned))
#  stats(secondDT.3.cleaned)
  #get.den <- density(secondDT.3.cleaned)
  #peak.quick(get.den$x, get.den$y)
 # length(secondDT.3.cleaned)



SP_count <-   num.of.DI.left - length(secondDT.3.cleaned)
}
SP <- list ("SP_mean" = SP_mean, "SP_std" = SP_std, "SP_count" = SP_count)
cleanedSample <- c(cleanedSample, SP)
#cleanedSample

##=====================================
# aneuploidy population of interest
##=====================================
aneup.pop <- secondDT.3.cleaned
plot(density(aneup.pop))
stats(aneup.pop)

aneu <- list ("AneuLeft" = aneup.pop)
cleanedSample <- c(cleanedSample, aneu)
cleanedSample

#cleanedSample$sample


##==========================
##  Saving the results
##==========================


storage.dir <- paste (root, "myGit/mixturemodel/cleanedData/OSCC/", sep = "")
file2save <- paste (storage.dir, "cleaned_", cleanedSample$sample, ".rda", sep="")
#file2save
save (cleanedSample, file = file2save)
}
