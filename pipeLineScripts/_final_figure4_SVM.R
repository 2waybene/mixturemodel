##=========================================================
##  File: _final_figure4_SVM.R
##  Author: Jianying Li
##  Comment: process the date new data only with two classes
##		 run through six different models and select
##		 SAM as a candidate
##		 The second parameter set seems work the best.
##  Extended (from _01): trying to predict "olk" class
##=========================================================


#===============================================
##	Set up os paths
##==============================================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

root <- windows
#root <- linux 
#root <- mac.os
##===============================================


library(caret)
##===================================
## set up working directory
## and getting the dataset
##===================================

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))


data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")

##	Trying to get the best figure 2, but have to use para2
##==============================================================
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")


#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")

##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL
 

  ##	Retain data ONLY with two classes




labels <- as.vector(data.2.classes$label)

data.2.classes <- data.2.classes[,-16]
data.2.classes <- cbind (data.2.classes, label=labels)
file2classes  <- data.2.classes



## create data partition
set.seed(123)
inTrainingSet <- createDataPartition(file2classes$label, p=.7, list=FALSE)
labelTrain <- file2classes[ inTrainingSet,]
labelTest <- file2classes[-inTrainingSet,]
nrow(labelTrain)
nrow(labelTest)


##=======================
##	rpart
##=======================
library(rpart)
library(partykit)


cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)


##====================================
##	SVM Example 
##====================================
set.seed(1)
svmTune1 <- train(label ~ ., 
		    data = labelTrain,

                 method = "svmRadial",
                 # The default grid of cost parameters go from 2^-2,
                 # 0.5 to 1, 
                 # We'll fit 9 values in that sequence via the tuneLength
                 # argument.
                 tuneLength = 9,
                 preProc = c("center", "scale"),
                 metric = "ROC",   
                 trControl = cvCtrl)





set.seed(1234)
svmTune4 <- train(label ~ ., 
		    data = labelTrain,

                 method = "svmRadial",
                 # The default grid of cost parameters go from 2^-2,
                 # 0.5 to 1, 
                 # We'll fit 9 values in that sequence via the tuneLength
                 # argument.
                 tuneLength = 9,
                 preProc = c("center", "scale"),
                 metric = "ROC",   
                 trControl = cvCtrl)






###############################################################
## Slide 119: Collecting Results With resamples
cvValues <- resamples(list(SVM1 = svmTune1, SVM4 = svmTune4))
summary(cvValues)
#trellis.par.set()
bwplot(cvValues, layout = c(3, 1))



