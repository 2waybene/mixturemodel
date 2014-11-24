
##=========================================================
##  File: _final_figure4.R
##  Author: Jianying Li
##  Comment: process the date new data only with two classes
##		 run through six different models and select
##		 SAM as a candidate
##		 The second parameter set seems work the best.
##  NOT current!! JYL Just for a record keeping
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


#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")

##	Trying to get the best figure 2, but have to use para2
##==============================================================
data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")


#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")

##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL
 

  ##	Retain data ONLY with two classes


  ##	Retain data ONLY with two classes

data.2.classes <- dataN0[-which (dataN0$label == "k"),]
data.k <- dataN0[which (dataN0$label == "k"),]


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


##=======================
##	rpart
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
##	SVM Example 
##====================================
#set.seed(1234)
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
##	partial logistic regression
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
##	randomForest
##========================

library(randomForest)
set.seed(1)
RFTune <- train(label ~.,data=labelTrain, 
				method = "rf",
		     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

##========================
##	knn
##========================
set.seed(1)
knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
 		     metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)

##========================
##	neural network
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

summary(cvValues)
#trellis.par.set()
bwplot(cvValues, layout = c(3, 1))






##==========================================================
##	Learning model training with Caret
##	http://topepo.github.io/caret/training.html#control
##==========================================================
trellis.par.set(caretTheme())


##================================
##	To reproduce teh figure 4
##================================
setwd(paste (root, "/myGit/mixturemodel/manuscript/_4ModlSelection/", sep=""))
getwd()

figure4 =  "FirstTrial_02.jpeg"

jpeg (figure4)
##=================
#	box plot
##=================
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))
dev.off()


=======
##=========================================================
##  File: _final_figure4.R
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
data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")

#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL
 

  ##	Retain data ONLY with two classes

data.2.classes <- dataN0[-which (dataN0$label == "k"),]
data.k <- dataN0[which (dataN0$label == "k"),]


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


##=======================
##	rpart
##=======================
library(rpart)
library(partykit)


cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)


#set.seed(1)
set.seed(12345)
rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)

##====================================
##	SVM Example 
##====================================
#set.seed(1)
set.seed(12345)
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
##	partial logistic regression
##========================================
set.seed(12345)
plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)


##========================
##	randomForest
##========================

library(randomForest)
set.seed(12345)
RFTune <- train(label ~.,data=labelTrain, 
				method = "rf",
		     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

##========================
##	knn
##========================
set.seed(12345)
knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
 		     metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)

##========================
##	neural network
##========================

library(MASS)
set.seed(12345)
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
set.seed(12345)
cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           plr = plrTune, nnet = nnetFit,
					knn = knnFit1, rrf = RFTune)
               )

###############################################################
## Slide 119: Collecting Results With resamples

summary(cvValues)
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))



##==========================================================
##	Learning model training with Caret
##	http://topepo.github.io/caret/training.html#control
##==========================================================
trellis.par.set(caretTheme())


##================================
##	To reproduce teh figure 4
##================================
setwd(paste (root, "/myGit/mixturemodel/manuscript/_4ModlSelection/", sep=""))
getwd()

figure4 =  "Figure4_final_11242014.jpeg"

jpeg (figure4)
##=================
#	box plot
##=================
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))
dev.off()


tiff(file="Figure4_final_11242014.tiff", height = 12, width = 17, units = 'cm',
     compression = "lzw", res = 300)
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))
dev.off()





>>>>>>> ecfaa69d1e40bc2ba4e111f8ff7a560d43067c22
