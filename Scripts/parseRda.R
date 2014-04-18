#  simulating DNA D.I. values
# Need to source ~/myGit/mixturemodel/Scripts/simDt_functions.R
#	normal population
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

dt.dir <- paste (root, "/myGit/mixturemodel/cleanedData/OSCC/", sep="")

files <- list.files (path = dt.dir, pattern=".rda")
files

aneuMax = 0;
aneuMin = 2.3;

for (i in 1:length(files))
{
   load(paste(dt.dir, files[i], sep=""))
   if (length(cleanedSample$AneuLeft) != 0)
  {
    if (aneuMax < max(cleanedSample$AneuLeft))
    {
       aneuMax = max(cleanedSample$AneuLeft)
    }
  }
}
r
aneuMax



