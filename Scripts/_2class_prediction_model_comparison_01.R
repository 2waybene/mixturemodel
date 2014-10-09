##=========================================================
##  File: _2class_prediction_model_comparison_01.R
##  Author: Jianying Li
##  Comment: process the date new data only with two classes
##		 run through six different models and select
##		 SAM as a candidate
##		 The second parameter set seems work the best.
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
#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")

#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

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
table(factor(data.2.classes$label))




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
rpartFull


library(partykit)
rpartFulla <- as.party(rpartFull)
plot(rpartFulla)


rpartPred <- predict(rpartFull, labelTest, type = "class")
confusionMatrix(rpartPred, labelTest$label)   # requires 2 factor vectors


##################################################################

cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
set.seed(1)
rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)

rpartTune


trellis.par.set(caretTheme())
plot(rpartTune, scales = list(x = list(log = 10)))



ggplot(rpartTune) +scale_x_log10()


###############################################################
## Slide 78: Test Set Results

rpartPred2 <- predict(rpartTune, labelTest)
confusionMatrix(rpartPred2, labelTest$label)


###############################################################
## Slide 79: Predicting Class Probabilities

rpartProbs <- predict(rpartTune, labelTest, type = "prob")
head(rpartProbs)



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

den.c <- density(svmPred.k.prob$c)
den.n <- density(svmPred.k.prob$n)

plot( den.c)
lines( den.n, col = "red")

file.olk[-which(svmPred.k.prob$c >0.5),]

plot(density(svmPred.k.prob$c))
plot(density(svmPred.k.prob$n))
pairs(svmPred.k.prob)



##======================================
##	partial logistic regression
##========================================
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

##========================
##	knn
##========================

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
                           plr = plrTune, nnet = nnetFit,
					knn = knnFit1, rrf = RFTune)
               )



###############################################################
## Slide 119: Collecting Results With resamples

summary(cvValues)
save(file = "two_class_comparison_6_model.Rdata", cvValues)

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



