##=========================================================
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

storage.dir <- paste (root, "myGit/mixturemodel/cleanedData/", sep = "")
cleanedFiles <- list.files (storage.dir, pattern = "rda")
cleanedFiles[1]

load (paste (storage.dir, cleanedFiles[1], sep=""))

str(cleanedSample)
