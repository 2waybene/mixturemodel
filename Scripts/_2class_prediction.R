<<<<<<< HEAD
##=========================================================
##  File: _2class_prediction.R
##  Author: Newly tested ONLY on normal and cancer sample
##  Author: Jianying
##=========================================================
library(caret)

file2classes  <- data.2.classes
file.olk <- data.k 


## create data partition

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
rpartFull <- rpart(label ~ ., data = labelTrain)
rpartFulla <- as.party(rpartFull)
plot(rpartFulla)

rpartPred <- predict(rpartFull, labelTest, type = "class")
confusionMatrix(rpartPred, labelTest$label)   # requires 2 factor vectors




cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
set.seed(1)
rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)
rpartPred2 <- predict(rpartTune, labelTest)
confusionMatrix(rpartPred2, labelTest$label)



##=======================
##	C5.0
##=======================
grid <- expand.grid(model = "tree",
                    trials = c(1:100),
                    winnow = FALSE)
set.seed(1)
c5Tune <- train(labelTrain, labelTrain$label,
                method = "C5.0",
                metric = "ROC",
                tuneGrid = grid,                    
                trControl = cvCtrl)

c5Pred <- predict(c5Tune, labelTest)
confusionMatrix(c5Pred, labelTest$label)


##====================================
##	SVM Example 
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

svmPred <- predict(svmTune, labelTest)
confusionMatrix(svmPred, labelTest$label)


##=============================
##	Test on olk sample
##=============================
dim(file.olk)
svmPred.k.prob  <- predict(svmTune, file.olk, type = "prob")
svmPred.k.prob 

file.olk[-which(svmPred.k.prob$c >0.5),]

plot(density(svmPred.k.prob$c))
plot(density(svmPred.k.prob$n))
pairs(svmPred.k.prob)

###############################################################
## Slide 117: A Few Other Models 

set.seed(1)
fdaTune <- train(label ~ ., data = labelTrain, 
                 method = "fda",
                 tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)


fdaPred <- predict(fdaTune, labelTest)
confusionMatrix(fdaPred, labelTest$label)

set.seed(1)
plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)

plrPred <- predict(plrTune, labelTest)
confusionMatrix(plrPred, labelTest$label)

##========================
##	randomForest
##========================

library(randomForest)

RFTune <- train(label ~.,data=labelTrain, 
				method = "rf",
		     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

RFPred2 <- predict(RFTune, labelTest)
confusionMatrix(RFPred2, labelTest$label)

RFtrain <- randomForest(label ~.,data=labelTrain )
RFPred <- predict(RFtrain, labelTest)
confusionMatrix(RFPred, labelTest$label)


rrfFit <- train(label ~ ., method = "RRF", data = labelTrain)
rrfPred <- predict(rrfFit, labelTest)
confusionMatrix(rrfPred, labelTest$label)



knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
 		     metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)


knnPred <- predict(knnFit1, labelTest)
confusionMatrix(knnPred , labelTest$label)


library(MASS)
nnetFit <- train(label ~.,data=labelTrain, 
                 method = "nnet",
                 preProcess = "range", 
                 tuneLength = 2,
 		     metric = "ROC",
                 trace = FALSE,
                 maxit = 100,
			trControl = cvCtrl)



nnetPred <- predict(nnetFit, labelTest)
confusionMatrix(nnetPred , labelTest$label)


###############################################################
## Slide 118: Collecting Results With resamples

cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           C5.0 = c5Tune, FDA = fdaTune, 
                           plr = plrTune, nnet = nnetFit,
					knn = knnFit1, rrf = RFTune)
               )



###############################################################
## Slide 119: Collecting Results With resamples

summary(cvValues)


##==========================================================
##	Learning model training with Caret
##	http://topepo.github.io/caret/training.html#control
##==========================================================
trellis.par.set(caretTheme())
plot(nnetFit)
plot(knnFit1)
plot(RFTune)
plot(svmTune)




##=================
#	box plot
##=================
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))


###############################################################
## Slide 120: Visualizing the Resamples

library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
splom(cvValues, metric = "ROC", pch = 16, 
      cex = .7, col = rgb(.2, .2, .2, .4), 
      pscales = 0)


###############################################################
## Slide 122: Visualizing the Resamples

trellis.par.set(caretTheme())
dotplot(cvValues, metric = "ROC")

###############################################################
## Slide 123: Comparing Models

rocDiffs <- diff(cvValues, metric = "ROC")
summary(rocDiffs)



###############################################################
## Slide 124: Visualizing the Differences

trellis.par.set(caretTheme())
dotplot(rocDiffs, metric = "ROC")



###############################################################
## Slide 128: Clustering the Models

plot(caret:::cluster.resamples(cvValues), sub = "", main = "")



##### BEGIN: train model - NaiveBayes >>>>>
##	FAILED
nbFit  <- train(
  label ~ .,
  data = labelTrain,
  method='nb',
  trControl=trainControl(method='cv',number=10)
  )

nbPred <- predict(nbFit , labelTest)
confusionMatrix(nbPred, labelTest$label)



=======
##=========================================================
##  File: _2class_prediction.R
##  Author: Newly tested ONLY on normal and cancer sample
##  Author: Jianying
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

data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL
 

  ##	Retain data ONLY with two classes

data.2.classes <- dataN0[-which (dataN0$label == "k"),]
data.k <- dataN0[which (dataN0$label == "k"),]



library(caret)

file2classes  <- data.2.classes
file.olk <- data.k 


## create data partition



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
rpartFull <- rpart(label ~ ., data = labelTrain)
rpartFulla <- as.party(rpartFull)
plot(rpartFulla)
confusionMatrix(rpartPred, labe
rpartPred <- predict(rpartFull, labelTest, type = "class")
lTest$label)   # requires 2 factor vectors




cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
set.seed(1)
rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)
rpartPred2 <- predict(rpartTune, labelTest)
confusionMatrix(rpartPred2, labelTest$label)



##=======================
##	C5.0
##=======================
grid <- expand.grid(model = "tree",
                    trials = c(1:100),
                    winnow = FALSE)
set.seed(1)
c5Tune <- train(labelTrain, labelTrain$label,
                method = "C5.0",
                metric = "ROC",
                tuneGrid = grid,                    
                trControl = cvCtrl)

c5Pred <- predict(c5Tune, labelTest)
confusionMatrix(c5Pred, labelTest$label)


##====================================
##	SVM Example 
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

svmPred <- predict(svmTune, labelTest)
confusionMatrix(svmPred, labelTest$label)


##=============================
##	Test on olk sample
##=============================
dim(file.olk)
svmPred.k.prob  <- predict(svmTune, file.olk, type = "prob")
svmPred.k.prob 

file.olk[-which(svmPred.k.prob$c >0.5),]

plot(density(svmPred.k.prob$c))
plot(density(svmPred.k.prob$n))
pairs(svmPred.k.prob)

###############################################################
## Slide 117: A Few Other Models 

set.seed(1)
fdaTune <- train(label ~ ., data = labelTrain, 
                 method = "fda",
                 tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)


fdaPred <- predict(fdaTune, labelTest)
confusionMatrix(fdaPred, labelTest$label)

set.seed(1)
plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)

plrPred <- predict(plrTune, labelTest)
confusionMatrix(plrPred, labelTest$label)

##========================
##	randomForest
##========================

library(randomForest)

RFTune <- train(label ~.,data=labelTrain, 
				method = "rf",
		     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

RFPred2 <- predict(RFTune, labelTest)
confusionMatrix(RFPred2, labelTest$label)

RFtrain <- randomForest(label ~.,data=labelTrain )
RFPred <- predict(RFtrain, labelTest)
confusionMatrix(RFPred, labelTest$label)


rrfFit <- train(label ~ ., method = "RRF", data = labelTrain)
rrfPred <- predict(rrfFit, labelTest)
confusionMatrix(rrfPred, labelTest$label)



knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
 		     metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)


knnPred <- predict(knnFit1, labelTest)
confusionMatrix(knnPred , labelTest$label)


library(MASS)
nnetFit <- train(label ~.,data=labelTrain, 
                 method = "nnet",
                 preProcess = "range", 
                 tuneLength = 2,
 		     metric = "ROC",
                 trace = FALSE,
                 maxit = 100,
			trControl = cvCtrl)



nnetPred <- predict(nnetFit, labelTest)
confusionMatrix(nnetPred , labelTest$label)


###############################################################
## Slide 118: Collecting Results With resamples

cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           C5.0 = c5Tune, FDA = fdaTune, 
                           plr = plrTune, nnet = nnetFit,
					knn = knnFit1, rrf = RFTune)
               )



###############################################################
## Slide 119: Collecting Results With resamples

summary(cvValues)


##==========================================================
##	Learning model training with Caret
##	http://topepo.github.io/caret/training.html#control
##==========================================================
trellis.par.set(caretTheme())
plot(nnetFit)
plot(knnFit1)
plot(RFTune)
plot(svmTune)




##=================
#	box plot
##=================
trellis.par.set()
bwplot(cvValues, layout = c(3, 1))


###############################################################
## Slide 120: Visualizing the Resamples

library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
splom(cvValues, metric = "ROC", pch = 16, 
      cex = .7, col = rgb(.2, .2, .2, .4), 
      pscales = 0)


###############################################################
## Slide 122: Visualizing the Resamples

trellis.par.set(caretTheme())
dotplot(cvValues, metric = "ROC")

###############################################################
## Slide 123: Comparing Models

rocDiffs <- diff(cvValues, metric = "ROC")
summary(rocDiffs)



###############################################################
## Slide 124: Visualizing the Differences

trellis.par.set(caretTheme())
dotplot(rocDiffs, metric = "ROC")



###############################################################
## Slide 128: Clustering the Models

plot(caret:::cluster.resamples(cvValues), sub = "", main = "")



##### BEGIN: train model - NaiveBayes >>>>>
##	FAILED
nbFit  <- train(
  label ~ .,
  data = labelTrain,
  method='nb',
  trControl=trainControl(method='cv',number=10)
  )

nbPred <- predict(nbFit , labelTest)
confusionMatrix(nbPred, labelTest$label)



>>>>>>> 0e2aeba706913c900210e72ae3d381a43cf26b6d
