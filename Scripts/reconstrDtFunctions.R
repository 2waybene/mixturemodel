  
##===================================================
# Main assumption in recontructing variables
##===================================================
para1 <- function()
{
  oneSampleRatio <- c(98,1.5,0.5)
  twoSampleRatio <- c(99.5, 0.5)
  triSampleRatio <- c(90,10)

  fakeSP_mean <- 2.0
  fakeSP_std  <- 0.3

  fake_aneu_mean <- 2.3
  fake_aneu_std  <- 0.3

  filler <- 0.00001
  paras <- list (oneSampleRatio = oneSampleRatio, 
                 twoSampleRatio = twoSampleRatio,
                 triSampleRatio = triSampleRatio,
                 fakeSP_mean = fakeSP_mean,
                 fakeSP_std  = fakeSP_std,
                 fake_aneu_mean = fake_aneu_mean,
                 fake_aneu_std  = fake_aneu_std,
                 filler = filler
  )
  return(paras)
}

para2 <- function()
{
  
  oneSampleRatio <- c(78,21.5,0.5)
  twoSampleRatio <- c(99.5, 0.5)
  triSampleRatio <- c(90,10)
  
  fakeSP_mean <- 2.0
  fakeSP_std  <- 0.3
  
  fake_aneu_mean <- 2.3
  fake_aneu_std  <- 0.3
  
  filler <- 0.00001
  paras <- list (oneSampleRatio = oneSampleRatio, 
                twoSampleRatio = twoSampleRatio,
                triSampleRatio = triSampleRatio,
                fakeSP_mean = fakeSP_mean,
                fakeSP_std  = fakeSP_std,
                fake_aneu_mean = fake_aneu_mean,
                fake_aneu_std  = fake_aneu_std,
                filler = filler
                )
  return(paras)
}


para3 <- function()
{
  
  oneSampleRatio <- c(99, 0.5, 0.5)
  twoSampleRatio <- c(99.5, 0.5)
  triSampleRatio <- c(90,10)
  
  fakeSP_mean <- 2.0
  fakeSP_std  <- 0.3
  
  fake_aneu_mean <- 2.3
  fake_aneu_std  <- 0.3
  
  filler <- 0.00001
  paras <- list (oneSampleRatio = oneSampleRatio, 
                 twoSampleRatio = twoSampleRatio,
                 triSampleRatio = triSampleRatio,
                 fakeSP_mean = fakeSP_mean,
                 fakeSP_std  = fakeSP_std,
                 fake_aneu_mean = fake_aneu_mean,
                 fake_aneu_std  = fake_aneu_std,
                 filler = filler
  )
  return(paras)
}

para4 <- function()
{
  
  oneSampleRatio <- c(99.8, 0.1, 0.1)
  twoSampleRatio <- c(99.8, 0.2)
  triSampleRatio <- c(90,10)
  
  fakeSP_mean <- 2.0
  fakeSP_std  <- 0.3
  
  fake_aneu_mean <- 2.3
  fake_aneu_std  <- 0.3
  
  filler <- 0.00001
  paras <- list (oneSampleRatio = oneSampleRatio, 
                 twoSampleRatio = twoSampleRatio,
                 triSampleRatio = triSampleRatio,
                 fakeSP_mean = fakeSP_mean,
                 fakeSP_std  = fakeSP_std,
                 fake_aneu_mean = fake_aneu_mean,
                 fake_aneu_std  = fake_aneu_std,
                 filler = filler
  )
  return(paras)
}


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

