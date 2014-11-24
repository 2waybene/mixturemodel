##=========================================================
##  File: Figure4.R
##  Author: Jianying Li
##  Initial tested as: _final_figure4.R
##  Comment: to produce figure 4
##=========================================================

#===============================================
##Set up os paths
##==============================================
#mac.os  <- "/Users/li11/"
#root <- mac.os

#linux   <- "~/"
#root <- linux 

windows <- "X:/"
root <- windows




##===============================================
library(caret)


##===================================
## set up working directory
## and getting the dataset
##===================================
setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))
data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")



##data cleaning
var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL


  ##Retain data ONLY with two classes
data.2.classes <- dataN0[-which (dataN0$label == "k"),]
dim(data.2.classes)
labels <- as.vector(data.2.classes$label)
data.2.classes <- data.2.classes[,-16]
data.2.classes <- cbind (data.2.classes, label=labels)
file2classes  <- data.2.classes




## create data partition
set.seed(1)
inTrainingSet <- createDataPartition(file2classes$label, p=.7, list=FALSE)
labelTrain <- file2classes[ inTrainingSet,]
labelTest <- file2classes[-inTrainingSet,]
nrow(labelTrain)
nrow(labelTest)


##=================================
##	Comparing models
##=================================

##=======================
##rpart
##=======================
library(rpart)
library(partykit)
cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
set.seed(1)
rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)
##====================================
##SVM Example 
##====================================
set.seed(1)
svmTune <- train(label ~ ., 
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


##======================================
##partial logistic regression
##========================================
set.seed(1)
plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)


##========================
##randomForest
##========================
library(randomForest)
set.seed(1)
RFTune <- train(label ~.,data=labelTrain, 
method = "rf",
     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)


##========================
##knn
##========================
set.seed(1)
knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
      	metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)


##========================
##neural network
##========================
library(MASS)
set.seed(1)
nnetFit <- train(label ~.,data=labelTrain, 
                 method = "nnet",
                 preProcess = "range", 
                 tuneLength = 2,
      metric = "ROC",
                 trace = FALSE,
                 maxit = 100,
trControl = cvCtrl)
###############################################################
## Slide 118: Collecting Results With resamples
set.seed(1)
cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           plr = plrTune, nnet = nnetFit,
				   knn = knnFit1, rrf = RFTune)
               )



###############################################################
## Slide 119: Collecting Results With resamples
setwd(paste (root, "/myGit/mixturemodel/manuscript/tables/", sep=""))
getwd()
sink ("Model_selection.txt")
summary(cvValues)
sink()

##===============================
#	For viewing purposes
##===============================

#trellis.par.set()
#bwplot(cvValues, layout = c(3, 1))


##==================

##================================
##To reproduce teh figure 4
##================================
setwd(paste (root, "/myGit/mixturemodel/manuscript/figures/", sep=""))
getwd()


tiff(file="Figure4_final_11242014.tiff", height = 12, width = 17, units = 'cm',
     compression = "lzw", res = 300)
bwplot(cvValues, layout = c(3, 1))
dev.off()



