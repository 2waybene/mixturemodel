##=============================================================
##	File name: R-peak-functions.R
##	Author   : Jianying Li
##	History  : Initial coded, 12/09/2013
##	Comment  : This custom function is a collection of peak
##		     finding functions. 
##=============================================================
source("x:/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")


###

#		peak.quick(x,y)
#				x 	- a numerical vector 
#				p(x)  - like normal denisty
#				return: peaks

#		peakfinder(x, peaks = 3)
#				x	- a numerical vector
#				peaks	- default find the top three peaks

#		findPeaks.mod ( x, thd =0)
#				Inherited from quantmod package with 1 off on the corrdinates
#				x	- a numerical vector
#				thd	- default find all peaks


##==================================##
##	Peak functions here		##
##==================================##


##	peak.quick	##

peak.quick <- function (x, y){

	return(x[which(diff(sign(diff(y)))==-2)])
}


peakfinder <- function(d, peaks = 3){
  dh <- hist(d,plot=FALSE)
  ins <- dh[["density"]]	#Newly update JYL
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins-(peaks -1),to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}


##	peak.dt.only	##

peak.dt.only <- function(x, halfWindowSize) {

  windowSize <- halfWindowSize * 2 + 1
  windows <- embed(x, windowSize)
  localMaxima <- max.col(windows, "first") == halfWindowSize + 1

  return(c(rep(FALSE, halfWindowSize), localMaxima, rep(FALSE, halfWindowSize)))
}


##	findPeaks in quantmod

library(quantmod)
findPeaks.mod <- function( x, threshold =0)
{
	p=findPeaks(x, threshold)
	return(p-1)
}



