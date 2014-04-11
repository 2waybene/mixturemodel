##======================================================
##	FileName: cleanFuncs.R
##  	clean the first population
##======================================================
cleanFirstPop <- function ( firstPeak, 
                            dt.first, 
                            dt.raw  )

{
##  Filter starts here...

first.den <- density(dt.first + firstPeak)
tempDen <-  first.den
tempDen$y <- first.den$y/sum(first.den$y)

##  Retain data on the right of the first peak
##			SAME AS
##	Remove data on the left  of the fist peak

	dt.right.of.peak <- dt.raw[which(dt.raw >=  peaks[1])]

##	Retain data on the right of max of the first population

	dt.right.of.peak.max <- dt.right.of.peak [which(dt.right.of.peak >=max(dt.first+peaks[1]))]

##	Data fall between the right of the first peak and the left of max of the first population
	dt.between.peak.max <- dt.right.of.peak [-which(dt.right.of.peak >=max(dt.first+peaks[1]))]

##	Now, remove the data according to the "estimated proportion"
##	Between two adjacent populations

	adjust = 0; 
	dt2filter <-dt.between.peak.max
#	str(dt2filter)
	for (i in 1:256)
	{
		temp = 0
	  	l.bound <- i + 255
  		h.bound <- i + 256
  		num.of.data <- ((tempDen$y[l.bound] + tempDen$y[h.bound])/2)*(length(dt.first))
  		candidate <- which(dt2filter > tempDen$x[l.bound] & dt2filter  < tempDen$x[h.bound])
   		if (length(candidate) >=1)
  		{
    			if (length(candidate) > floor(num.of.data))
    			{
      			temp = num.of.data - floor(num.of.data)
      			data2exclude <- sample(candidate, floor(num.of.data))
      			if (length(data2exclude) >=1 )
      			{
        				dt2filter <- dt2filter[-data2exclude]	
        				adjust = adjust + temp
      			}
    			}else{
      			dt2filter <- dt2filter[-candidate]
    			}
  		}
	}
	num2salvage <- sample (c(1:length(dt2filter)), ceiling(adjust)) 	#FIXME: Manully fixting
	dt2filter <- dt2filter[-num2salvage]
	dt.retain <- c(dt.right.of.peak.max, dt2filter)
	return (dt.retain)

}