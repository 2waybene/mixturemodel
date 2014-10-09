##============================================
##  File: predModHongMod_v1.R
##  Author: Hong Xu
##  Modified from preModHong_02.R by Hong Xu
##  then modified to predModHongMod.R
##  Now, it becomes predModHongMod_v1.R
##============================================
library(caret)
library(pROC)
library(Metrics)


#===========================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- windows
root <- mac.os
#==========================






setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))

##	param1
#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
#setwd(paste (root, "/myGit/mixturemodel/modeling/", sep=""))
#sink ("log_param1.txt")

##	param2
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")
#setwd(paste (root, "/myGit/mixturemodel/modeling/", sep=""))
#sink ("log_param2.txt")


##	param3
data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")
setwd(paste (root, "/myGit/mixturemodel/modeling/", sep=""))
#sink ("log_param3.txt")


##	param4
#data <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")
#setwd(paste (root, "/myGit/mixturemodel/modeling/", sep=""))
#sink ("log_param4.txt")                   
op <- par()
par (mfrow = c(1,3))
data.c <- data[which(data$labe =="c"),]
boxplot(data.c[,c(-1, -17, -18)], main = "ExGCRn -- OSCC")



data.n <- data[which(data$labe =="n"),]
boxplot(data.n[,c(-1, -17, -18)], main = "ExGCRn -- Normal")


data.k <- data[which(data$labe =="k"),]
boxplot(data.k[,c(-1, -17, -18)], main = "ExGCRn -- OLK")
par(op)
##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL


