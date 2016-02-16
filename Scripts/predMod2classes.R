##=========================================================
##  File: predMod2classes.R
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
table(factor(data.2.classes$label))




#levels(data.2.classes$label) <- factor (data.2.classes$label)

                  

## create data partition

inTrainingSet <- createDataPartition(data.2.classes$label, p=.7, list=FALSE)
labelTrain <- data.2.classes[ inTrainingSet,]
labelTest <- data.2.classes[-inTrainingSet,]
nrow(labelTrain)
nrow(labelTest)






##### END: data partition <<<<<
##### BEGIN: tune the parameters >>>>>



library(rpart)
rpart1 <- rpart(label ~ ., data = labelTrain, 
                control = rpart.control(maxdepth = 2))
rpart1

###############################################################
## Slide 58: Visualizing the Tree

library(partykit)
rpart1a <- as.party(rpart1)
plot(rpart1a)


###############################################################
## Slide 61: The Final Tree

rpartFull <- rpart(label ~ ., data = labelTrain)

###############################################################
## Slide 62: The Final Tree

rpartFull
#plot(rpartFull)
###############################################################
## Slide 64: The Final rpart Tree

rpartFulla <- as.party(rpartFull)
plot(rpartFulla)


###############################################################
## Slide 65: Test Set Results

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


###############################################################
## Slide 73: train Results

rpartTune


###############################################################
## Slide 75: Resampled ROC Profile

trellis.par.set(caretTheme())
plot(rpartTune, scales = list(x = list(log = 10)))

###############################################################
## Slide 76: Resampled ROC Profile

ggplot(rpartTune) +scale_x_log10()


###############################################################
## Slide 78: Test Set Results

rpartPred2 <- predict(rpartTune, labelTest)
confusionMatrix(rpartPred2, labelTest$label)


###############################################################
## Slide 79: Predicting Class Probabilities

rpartProbs <- predict(rpartTune, labelTest, type = "prob")
head(rpartProbs)


###############################################################
## Slide 80: Creating the ROC Curve

library(pROC)
rpartROC <- roc(labelTest$label, rpartProbs[, "c"], 
                levels = rev(levels(labelTest$label)))

plot(rpartROC, type = "S", print.thres = .5)


##	rpart done here



## Slide 91: Tuning the C5.0 Model

grid <- expand.grid(model = "tree",
                    trials = c(1:100),
                    winnow = FALSE)
set.seed(1)
c5Tune <- train(labelTrain, labelTrain$label,
                method = "C5.0",
                metric = "ROC",
                tuneGrid = grid,                    
                trControl = cvCtrl)



###############################################################
## Slide 92: Model Output

c5Tune

###############################################################
## Slide 93: Boosted Tree Resampling Profile

ggplot(c5Tune)

###############################################################
## Slide 94: Test Set Results

c5Pred <- predict(c5Tune, labelTest)
confusionMatrix(c5Pred, labelTest$label)

#dim(data.k)
#lab <- c (rep("c", 82))
#levels(as.factor(lab)) <- levels(as.factor(c("c", "n")))
#levels(data.k$label)


#c5Pred <- predict(c5Tune, data.k)
#confusionMatrix(c5Pred, data.k$label)

###############################################################
## Slide 95: Test Set ROC Curve

c5Probs <- predict(c5Tune, labelTest, type = "prob")
head(c5Probs, 3)

library(pROC)
c5ROC <- roc(predictor = c5Probs$n,
             response = labelTest$label,
             levels = rev(levels(labelTest$label)))


###############################################################
## Slide 96: Test Set ROC Curve

c5ROC

plot(rpartROC, type = "S")
plot(c5ROC, add = TRUE, col = "#9E0142")

###############################################################
## Slide 98: Test Set Probabilities

histogram(~c5Probs$n|labelTest$label, xlab = "Probability of Poor Segmentation")

###############################################################
## Slide 109: SVM Example 

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


###############################################################
## Slide 110: SVM Example

svmTune


###############################################################
## Slide 111: SVM Example

svmTune$finalModel



###############################################################
## Slide 112: SVM ROC Profile

plot(svmTune, metric = "ROC", scales = list(x = list(log = 2)))


###############################################################
## Slide 113: Test Set Results

svmPred <- predict(svmTune, labelTest)
confusionMatrix(svmPred, labelTest$label)

###############################################################
## Slide 114: Test Set ROC Curves

svmROC <- roc(labelTest$label, 
              predict(svmTune, labelTest, type = "prob")[, "n"],
              levels = rev(levels(labelTest$label)))
plot(rpartROC, type = "S")
plot(c5ROC, add = TRUE, col = "#9E0142")
plot(svmROC, add = TRUE, col = "#3288BD")

legend(.5, .5, c("CART", "C5.0", "SVM"), 
       col = c("black", "#9E0142", "#3288BD"), 
       lty = rep(1, 3))



###############################################################
## Slide 117: A Few Other Models 

set.seed(1)
fdaTune <- train(label ~ ., data = labelTrain, 
                 method = "fda",
                 tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

set.seed(1)
plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)


###############################################################
## Slide 118: Collecting Results With resamples

cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           C5.0 = c5Tune, FDA = fdaTune, 
                           logistic = plrTune))



###############################################################
## Slide 119: Collecting Results With resamples

summary(cvValues)



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




#################################################################


library(randomForest)
modelRF <- randomForest(label ~.,data=labelTrain )



## control:
# resampling technique: 5-repeat 10-fold cross-validation
# performance metrics: ROC AUC curve

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE)


##### END: tune the parameters <<<<<
##### BEGIN: train model - svm >>>>>
set.seed(1024)
svmFit <- train(label ~ ., data = labelTrain,
                ## training model: svm >>>
                method = "svmRadial",
                metric = "ROC",
                tuneLength = 10,
                trControl = ctrl)
## prediction
svmPredProb <- predict(svmFit, labelTest ,type = "prob")

svmPred <- predict(svmFit, labelTest )
str(svmPred)

confusionMatrix(svmPred, labelTest$label)



## prediction on olk!!


svmProbs <- predict(svmFit, data.k , type = "prob")

str(svmProbs)

hist(svmProbs$c)

confusionMatrix(svmPred, labelTest$label)


## Confusion Matrix and Statistics
##
## Reference
## Prediction c k n
## c 8 2 0
## k 0 2 0
## n 2 2 30
##
## Overall Statistics
##
## Accuracy : 0.87
## 95% CI : (0.737, 0.951)
## No Information Rate : 0.652
## P-Value [Acc > NIR] : 0.000839
##
## Kappa : 0.72
## Mcnemar's Test P-Value : 0.111610
##
## Statistics by Class:
##
## Class: c Class: k Class: n
## Sensitivity 0.800 0.3333 1.000
## Specificity 0.944 1.0000 0.750
## Pos Pred Value 0.800 1.0000 0.882
## Neg Pred Value 0.944 0.9091 1.000
## Prevalence 0.217 0.1304 0.652
## Detection Rate 0.174 0.0435 0.652
## Detection Prevalence 0.217 0.0435 0.739
## Balanced Accuracy 0.872 0.6667 0.875
## plot
for(stat in c( 'Accuracy',
                'Kappa',
                'AccuracyLower',
                'AccuracyUpper',
  'AccuracyPValue',
  'Sensitivity',
  'Specificity',
  'Pos_Pred_Value',
  'Neg_Pred_Value',
  'Detection_Rate',
  'ROC',
  'logLoss')) 
{
  print(plot(svmFit, metric=stat))
}
##### END: train model - svm <<<<<
##### BEGIN: train model - rf >>>>>
rfFit <- train(label ~ ., method = "rf", data = labelTrain)
rfPred <- predict(rfFit, labelTest)
confusionMatrix(rfPred, labelTest$label)



## Confusion Matrix and Statistics
##
## Reference
## Prediction c k n
## c 8 0 0
## k 0 3 0
## n 2 3 30
##
## Overall Statistics
##
## Accuracy : 0.891
## 95% CI : (0.764, 0.964)
## No Information Rate : 0.652
## P-Value [Acc > NIR] : 0.000217
##
## Kappa : 0.762
## Mcnemar's Test P-Value : NA
##
## Statistics by Class:
##
## Class: c Class: k Class: n
## Sensitivity 0.800 0.5000 1.000
## Specificity 1.000 1.0000 0.688
## Pos Pred Value 1.000 1.0000 0.857
## Neg Pred Value 0.947 0.9302 1.000
## Prevalence 0.217 0.1304 0.652
## Detection Rate 0.174 0.0652 0.652
## Detection Prevalence 0.174 0.0652 0.761
## Balanced Accuracy 0.900 0.7500 0.844
##### END: train model - rf <<<<<
##### BEGIN: train model - rf >>>>>
rrfFit <- train(label ~ ., method = "RRF", data = labelTrain)
rrfPred <- predict(rrfFit, labelTest)
confusionMatrix(rrfPred, labelTest$label)
## Confusion Matrix and Statistics
##
## Reference
## Prediction c k n
## c 8 0 0
## k 0 4 1
## n 2 2 29
##
## Overall Statistics
##
## Accuracy : 0.891
## 95% CI : (0.764, 0.964)
## No Information Rate : 0.652
## P-Value [Acc > NIR] : 0.000217
##
## Kappa : 0.774
## Mcnemar's Test P-Value : NA
