#  simulating DNA D.I. values
# Need to source ~/myGit/mixturemodel/Scripts/simDt_functions.R
#	normal population
##======================================================
#	File name: recontrData-v1.R
#	Author: Jianying Li
#	Comment: used to convert OSCC to 
##=====================================================
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

source (paste (root, "myGit/mixturemodel/Scripts/cleaningFuncs.R", sep = ""))
source (paste (root, "myGit/mixturemodel/Scripts/simDt_functions.R", sep = ""))

dt.dir <- paste (root, "/myGit/mixturemodel/cleanedData/Normal/", sep="")

files <- list.files (path = dt.dir, pattern=".rda")
files

aneuMax = 0;
aneuMin = 2.3;

reconStrDt <- list(c())

dt.return <- ""

for (k in 1:length(files))
{
 k =1
#i =10

   load(paste(dt.dir, files[k], sep=""))

str(cleanedSample)
str(cleanedSample$AneuLeft)

	#Make it 8 if greater than 8
	if (cleanedSample$AneuLeft == "")
	{
		popNum = 2
	}
	
	if (cleanedSample$SP_count == "")
	{
		popNum = 1
	}

	cleanedSample$AneuLeft[which(cleanedSample$AneuLeft > 8)] <- 8 

	#if (is.na(cleanedSample$SP_count1))

	ratio <- cleanedSample$FP_count/cleanedSample$SP_count

	w.norm <- ratio/(1+ratio)
	w.mito <- 1/(1+ratio)
	x <- seq(0,2.3, by=(2.3/512))
	x <- x[-1]
	length(x)
	y1 <- w.norm*P(x, cleanedSample$FP_mean, cleanedSample$FP_std)
	y2 <- w.mito*P(x, cleanedSample$SP_mean, cleanedSample$SP_std)

	y = y1 + y2
length(y)
length(x)


plot(x,y/sum(y), type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")


	prob.y <- c()
	pdf.y <- y/sum(y)
	prob.y[1] <- pdf.y[1]

	for (i in 2:length(y))
	{
		temp.prob <- pdf.y[i]
		prob.y[i] <- prob.y[i-1] + temp.prob
	}


	if (length(cleanedSample$AneuLeft) < 5)
	{
		num2Recontr <- 1000
	}else{
		num2Recontr <- 9*numOfAneu
	}

	simDt <- c()
	seed = 12345
	for (i in 1:num2Recontr)
	{
		x1 <- runif(1, 0, 1)
		index <- which(prob.y < (x1 + 0.002) & prob.y > (x1 - 0.002))
		if (length(index) < 1)
		{
			index <- which(prob.y < (x1 + 0.02) & prob.y > (x1 - 0.02))
		}
		if (length(index) < 1)
		{
			index <- which(prob.y < (x1 + 0.2) & prob.y > (x1 - 0.2))
		}
		temp <- sample( index, 1)
		simDt[i] <- x[temp]
	}


stats(simDt)

	if (popNum ==3)
	{
		simDt <- c(simDt , cleanedSample$AneuLeft)
	}


plot(density(simDt))

	bk = floor(max(simDt))*2
	den <- (hist(simDt, breaks=bk)$density)
	for (m in length(den):16)
	{
		den[m] <- 0.00001
	}

#	reconStrDt[[k]]<- list (name=cleanedSample$sample, dat = den)

	dt.temp <- as.data.frame(den)
	colnames(dt.temp) <- cleanedSample$sample
	dt.return <- cbind(dt.return, dt.temp)
	
}

normal.temp <- dt.return[,-1]


dim(t(normal.temp))[1]
label <- rep("n", dim(t(normal.temp))[1])
normal.out <- cbind(t(as.data.frame(normal.temp)), as.data.frame(label))


