##====================================================
##	File:   simulate_DI.R
##	Author: Jianying Li
##	Comment: Learning mixture model for Aneuploidy
##		   study for OSCC
##====================================================


##	Technical hurdles
##	Unbalanced data distribution ??



source("x:/myGit/mixturemodel/ScriptsDraft/mixtureModelFunctions.R")
#source("/Users/li11/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")



##=====Sample data
macFileDir <- "/Users/li11/myGit/mixturemodel/data/"
winFileDir <- "x:/myGit/mixturemodel/data/"

f <- "fmd_DT_raw.txt"
f <- "oscc-olk1_parsed.txt"
f_IN <- paste (file.dir, f , sep="")
f_IN <- paste (macFileDir, f , sep="")
f_IN <- paste (winFileDir, f , sep="")
##=====End of settings

dt <- read.table(f_IN, header= F, sep = "\t")

str(dt)
temp.cv <- sd(dt$V1)/mean(dt$V1)
temp.cv

##=====End

##======================================================================
##  Based on our prior knowledge, we can focus on 
##  three clusters of cell populations:
##  Cluster one, with mean DI = 1.001
##  Cluster two, with mean DI = 2.001
##  Cluster three, with minimum DI > 2.300, assigning 2.300 as mean
##  Now, let's assume that C.V. is same 
##======================================================================

mean <- c(1.001, 2.002, 2.300)
#mean <- c(1.001, 2.002, 3.600)
#mean <- c(0.936, 1.692, 3.600)
sigma <- c(0.19, 0.25,0.5) 
sigma <- mean*temp.cv

weight <- c(0.893, 0.092, 0.05)
#weight <- c(0.425, 0.425, 0.05)


Delta <- 0.01
x <- seq(0,7, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3


#integrate(P, 0,7, mean[1], sigma[1])
#integrate(P, 0,7, mean[2], sigma[2])
#integrate(P, 0,7, mean[3], sigma[3])

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



