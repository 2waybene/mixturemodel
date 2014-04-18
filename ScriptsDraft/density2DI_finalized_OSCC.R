library(Rlab)

##  OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

##=============================================
# Three criteria
# aneuploidy: > 2.3c --> three populations
# mitotic: >1.7c --> two populations
# normal??

aneuThresh = 2.3
mitoThresh = 1.7
##=============================================
##  Read in data
##=============================================
#root <- windows
root <- mac.os

source (paste (root, "myGit/mixturemodel/Scripts/cleaningFuncs.R", sep = ""))
rawFiles <- list.files (paste(root, "myGit/mixturemodel/data/dt_01232014/OSCC/",sep=""), pattern = "csv")
rawFiles

i = 3
i = 4 

for (i in 1:length(rawFiles))
{
  i=6
print (rawFiles[i])
#cat ("\n")
#}

fileName <- paste("myGit/mixturemodel/data/dt_01232014/OSCC/", rawFiles[i], sep ="")


##========================================================
#  On a real D.I. value data
##=========================================================



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

##===================================================================
get.den <- density(as.vector(dt$DNA_Index))
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
dt.raw <- as.vector (dt$DNA_Index)
firstDT <- getPopWIndex (dt.raw, index)

##  Save first population dt
FP_dt_primary <- firstDT + peaks[index]

#FP_mean <- mean(firstDT + peaks[index])
#FP_std <- sd(firstDT + peaks[index])

dt.cleaned <- cleanFirstPop(peaks[index], firstDT, dt.raw)
#plot(density(dt.cleaned))

#str(dt.cleaned)

##================================
##  Extra cleaning if necessary
##================================
firstDT <- getFirstPop(dt.cleaned)
#plot(density(firstDT))
peaks <- peak.quick(density(dt.cleaned)$x, density(dt.cleaned)$y)
peaks

##=========================================
##  Follow the same protocol, but just 
##  carry out one more cleaning cycle
##  if there is any peaks less than 1.2
##===========================================

##Need to add the "cleaned back to population one"!!

if (peaks[1] < 1.2)
{
  dt.another.clean <- cleanFirstPop(peaks[1], firstDT, dt.cleaned)
#  plot(density(dt.another.clean))
  dt.1pop.cleaned <- dt.another.clean
}else{
  dt.1pop.cleaned <- dt.cleaned
}

FP_mean <- mean(FP_dt_primary)
FP_std <- sd(FP_dt_primary)
  FP_count <- (length(dt.raw) - length(dt.1pop.cleaned))

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
