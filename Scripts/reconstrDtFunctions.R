
##===================================================
# Main assumption in recontructing variables
##===================================================
oneSampleRatio <- c(98,1.5,0.5)
twoSampleRatio <- c(99.5, 0.5)
triSampleRatio <- c(90,10)

fakeSP_mean <- 2.0
fakeSP_std  <- 0.3

fake_aneu_mean <- 2.3
fake_aneu_std  <- 0.3

filler <- 0.00001

##===================================
#  Test functions here
##===================================


##===================================
#	Functions here
##===================================

getSimNum <- function (num = 100, low = 2, high = 8)
{
  mean.norm <- c()
  for (i in 1:num)
  {
    mean.norm[i] <- runif (1, low, high)
  }
  return (mean.norm)
}

getMeCDF <- function (x, y)
{
  range <- range(x)	
  bin <- (max(x) - min(x))/50
  for (i in (1:50))
  {
    
  }
  
  
}

