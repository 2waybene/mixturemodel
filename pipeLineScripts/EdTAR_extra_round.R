##========================
##	Another round??
##========================

dt <- dt.cleaned



## determine how many families are we dealing with
	numOfFamily <-  1 # minimun one family

	if (length(which(as.vector(dt) > aneuThresh)) >= 1)
	{
  		numOfFamily = 3
	}else if (length(which(as.vector(dt) > mitoThresh)) > 1)
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

	dt.raw <- as.vector (dt)
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
	if (peaks[length(which(peaks<1)) + 1] < normMax ) 
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

plotDensity(as.data.frame(dt.cleaned), xlab = "DNA Index value")

	get.den <- density(dt.cleaned)
	peaks <- peak.quick (get.den$x, get.den$y)
peaks

