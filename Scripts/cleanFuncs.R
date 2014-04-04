##======================================================
##  clean the first population
##======================================================
cleanFirstPop <- function ( firstPeak, 
                            dt.first, 
                            dt.raw  )

{
firstPeak <- peaks[1]
dt.first <- dt.first
dt.raw <- dt.raw


first.den <- density(dt.first + firstPeak)
tempDen <- first.den
tempDen$y <- first.den$y/sum(first.den$y)


# Data on the right side of first peak

dt.peak.right <- dt.raw[which(dt.raw >=  firstPeak)]
dt.no.first.pop  <- dt.peak.right [which(dt.peak.right >= max(dt.first+firstPeak))]
str(dt.no.first.pop)

##	Data fall between the right of the first peak and the left of max of the first population
dt.2.work <- dt.peak.right[-which(dt.peak.right >=max(dt.first+ firstPeak))]
summary(dt.raw.flt.2)
str(dt.raw.flt.2)


##	Now, remove the data according to the "estimated proportion"
##	Between two adjacent populations


adjust = 0; 
dt2filter <- dt.raw.flt.2
str(dt2filter)
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

str(dt2filter)
adjust

num2salvage <- sample (c(1:length(dt2filter)), ceiling(adjust)) 	#FIXME: Manully fixting
dt2filter <- dt2filter[-num2salvage]
#num2salvage <- sample (c(1:length(dt2filter)), 138)			#FIXME: Manully fixting!!!
#dt2filter <- dt2filter[-num2salvage]
str(dt2filter)

dt.raw.02 <- c(dt.raw.02, dt2filter)
str(dt.raw.02)
plot(density(dt.raw.02), main ="Density after removing the first population")

dt.second.pop <- dt.raw.02
str(dt.second.pop)


}