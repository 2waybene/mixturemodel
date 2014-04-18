##=============================================================
##	File name: R-peak-functions-test.R
##	Author   : Jianying Li
##	History  : Initial coded, 12/09/2013
##	Comment  : It is important to compare different
##		     peak functions
##=============================================================
source("x:/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")

##	Three (really two) peaks data
mean <- c(1.001, 2.002, 2.300)
sigma <- c(0.19, 0.25,0.5) 

weight <- c(0.893, 0.092, 0.05)
#weight <- c(0.425, 0.425, 0.05)

Delta <- 0.01
x <- seq(0,7, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3

plot(x,y) 	## Only TWO visible peaks

##	A better plot

par(mfrow=c(1,1))
plot(x,y, type="l", lwd=3,
     main="Simulated D.I. values: three-category groups",
     xlab="D.I. value", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")

lgd = c("Mixture", "Normal, mean = 1.001", "Mitotic,  mean = 2.002", "Aneuploid, mean = 2.300")
legend ("topright", lgd, text.col = c("black", "red", "green", "blue"))


##	Calling peaks...

peakx <- peak.quick(x,y)
peakx

peakx <- peak.quick(x,y1)
peakx #found the obvious peaks, but miss the trivil one(s)!!



p=findPeaks(x, 1)


peakx <- peak.dt.only(x,0.1)

which(peakx == TRUE)


##==================================##
##	Peak functions here		##
##==================================##


##	peak.quick	##

peak.quick <- function (x, y){

	return(x[which(diff(sign(diff(y)))==-2)])
}

##	peak.dt.only	##

peak.dt.only <- function(x, halfWindowSize) {

  windowSize <- halfWindowSize * 2 + 1
  windows <- embed(x, windowSize)
  localMaxima <- max.col(windows, "first") == halfWindowSize + 1

  return(c(rep(FALSE, halfWindowSize), localMaxima, rep(FALSE, halfWindowSize)))
}

x <- c(1,3,1,3,1)

peaks(x, 1)


data <- c(rnorm(100,mean=20),rnorm(100,mean=12))

peakfinder <- function(d){
  dh <- hist(d,plot=FALSE)
  ins <- dh[["density"]]	#Newly update JYL
#  ins <- dh[["intensities"]] #Newly update JYL
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins-2,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

peaks <- peakfinder(data)
hist(data)
sapply(peaks,function(x) abline(v=x,col="red"))



peaks <- peaks(x)
plot(x, type = "l")
points(which(peaks),x[which(peaks)])

peaks (data,1)
length(density(data))
density(data)

hist(data)
sapply(peaks,function(x) abline(v=x,col="red"))


aa=100:1
bb=sin(aa/3)
cc=aa*bb
plot(cc, type="l")
p=findPeaks.mod(cc)
points(p, cc[p])
p

