library(Rlab)

##  OS specific directories:

mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"


##=============================================
##  Read in data
##=============================================
#root <- windows
root <- mac.os

rawFiles <- list.files (paste(root, "myGit/mixturemodel/data/dt_01232014/",sep=""), pattern = "csv")
rawFiles

i =2 
rawFiles[i]

fileName <- paste("myGit/mixturemodel/data/dt_01232014/", rawFiles[i], sep ="")

##========================================================
#  On a real D.I. value data
##=========================================================
source (paste (root, "myGit/mixturemodel/Scripts/cleaningFuncs.R", sep = ""))


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
get.den <- density(as.vector(dt$DNA_Index))
plot(get.den)
peaks <- peak.quick (get.den$x, get.den$y)
peaks


dt.raw <- as.vector (dt$DNA_Index)
summary(dt.raw)
firstDT <- getFirstPop(dt.raw)

##  Save first population stats
FP_mean <- mean(firstDT + peak1)
FP_std <- sd(firstDT + peak1)

#FP_count <- length(firstDT + peak1)
#FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
#cleanedSample <- c(cleanedSample, FP)
#cleanedSample
peak1
str(dt.raw)
str(firstDT)
dt.cleaned <- cleanFirstPop(peak1, firstDT, dt.raw)
plot(density(dt.cleaned))

str(dt.cleaned)

##================================
##  Extra cleaning if necessary
##================================
firstDT <- getFirstPop(dt.cleaned)
plot(density(firstDT))
peaks <- peak.quick(density(dt.cleaned)$x, density(dt.cleaned)$y)

dt.another.clean <- cleanFirstPop(peaks[1], firstDT, dt.cleaned)
plot(density(dt.another.clean))

dt.1pop.cleaned <- dt.another.clean

FP_count <- (length(dt.raw) - length(dt.1pop.cleaned))

#FP_count <- length(firstDT + peak1)
FP <- list ("FP_mean" = FP_mean, "FP_std" = FP_std, "FP_count" = FP_count)
cleanedSample <- c(cleanedSample, FP)
cleanedSample


##================================
##  Get second population stats
##================================

num.of.DI.left <- length(dt.1pop.cleaned)

get.den <- density(dt.1pop.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
peaks


secondDT <- getSecondPop(dt.1pop.cleaned)
plot(density(secondDT))
##  Save first population stats
SP_mean <- mean(secondDT + peaks[2])
SP_std <- sd(secondDT + peaks[2])
#FP_count <- length(secondDT + peaks[2])
#FP <- list ("SP_mean" = FP_mean, "SP_std" = FP_std, "SP_count" = FP_count)
#cleanedSample <- c(cleanedSample, FP)
#cleanedSample

plot(density(secondDT + peaks[2]))


secondDT.cleaned <- cleanFirstPop(peaks[2],  secondDT, dt.1pop.cleaned)
str(secondDT.cleaned)
plot(density(secondDT.cleaned))


##==================================
##  Need another round of cleaning
##==================================

get.den <- density(secondDT.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
peaks

stats (secondDT.cleaned)

secondDT.1 <- getFirstPop(secondDT.cleaned)
plot(density(secondDT.1))
##  Save first population stats
#FP_mean <- mean(secondDT + peaks[2])
#FP_std <- sd(secondDT + peaks[2])
#FP_count <- length(secondDT + peaks[2])
#FP <- list ("SP_mean" = FP_mean, "SP_std" = FP_std, "SP_count" = FP_count)
#cleanedSample <- c(cleanedSample, FP)
#cleanedSample

plot(density(secondDT.1 + peaks[1]))


secondDT.2.cleaned <- cleanFirstPop(peaks[1],  secondDT.1, secondDT.cleaned)
str(secondDT.2.cleaned)
plot(density(secondDT.2.cleaned))
stats (secondDT.2.cleaned)


##  third round??

get.den <- density(secondDT.2.cleaned)
peaks <- peak.quick (get.den$x, get.den$y)
peaks

secondDT.2.sub <- getFirstPop(secondDT.2.cleaned)
secondDT.3.cleaned <- cleanFirstPop(peaks[1], secondDT.2.sub , secondDT.2.cleaned)
plot(density(secondDT.3.cleaned))
stats(secondDT.3.cleaned)
get.den <- density(secondDT.3.cleaned)
peak.quick(get.den$x, get.den$y)
length(secondDT.3.cleaned)



SP_count <-   num.of.DI.left - length(secondDT.3.cleaned)
SP <- list ("SP_mean" = SP_mean, "SP_std" = SP_std, "SP_count" = SP_count)
cleanedSample <- c(cleanedSample, SP)
cleanedSample

##=====================================
# aneuploidy population of interest
##=====================================
aneup.pop <- secondDT.3.cleaned
plot(density(aneup.pop))
stats(aneup.pop)

aneu <- list ("AneuLeft" = aneup.pop)
cleanedSample <- c(cleanedSample, aneu)
cleanedSample

cleanedSample$sample


##==========================
##  Saving the results
##==========================


storage.dir <- paste (root, "myGit/mixturemodel/cleanedData/", sep = "")
file2save <- paste (storage.dir, "cleaned_", cleanedSample$sample, ".rda", sep="")
file2save
save (cleanedSample, file = file2save)

