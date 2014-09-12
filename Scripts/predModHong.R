##  File: predModHong.R
##  Author: Hong Xu
library(caret)
library(pROC)
library(ROCR)
library(Metrics)

#====================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- windows
root <- mac.os


##### REF: http://stats.stackexchange.com/questions/31579/what-is-the-optimal-k-for-the-k-nearest-neighbour-classifier-on-the-iris-dat
# https://gist.github.com/zachmayer/3061272
#Multi-Class Summary Function
#Based on caret:::twoClassSummary

require(compiler)
multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL)
{
                            #Load Libraries
    require(Metrics)
    require(caret)
    #Check data
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
      stop("levels of observed and predicted data do not match")
      #Calculate custom one-vs-all stats for each class
     prob_stats <- lapply(levels(data[, "pred"]), function(class)
       {
          #Grab one-vs-all data for the class
           pred <- ifelse(data[, "pred"] == class, 1, 0)
           obs <- ifelse(data[, "obs"] == class, 1, 0)
           prob <- data[,class]
         #Calculate one-vs-all AUC and logLoss and return
           cap_prob <- pmin(pmax(prob, .000001), .999999)
           prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
           names(prob_stats) <- c("ROC", "logLoss")
           return(prob_stats)
      })
      prob_stats <- do.call(rbind, prob_stats)
      rownames(prob_stats) <- paste( "Class:" , levels(data[, "pred"]))
    
      #Calculate confusion matrix-based statistics
      CM <- confusionMatrix(data[, "pred"], data[, "obs"])
                            #Aggregate and average class-wise stats
                            #Todo: add weights
      class_stats <- cbind(CM$byClass, prob_stats)
      class_stats <- colMeans(class_stats)
                            #Aggregate overall stats
      overall_stats <- c(CM$overall)
    
      #Combine overall with class-wise stats and remove some stats we don't want
      stats <- c(overall_stats, class_stats)
      stats <- stats[! names(stats) %in% c("AccuracyNull","Prevalence", "Detection Prevalence")]
                            
    #Clean names and return
                  
    names(stats) <- gsub('[[:blank:]] +', '_' , names(stats))
    return(stats)
})

## Note: no visible binding for global variable 'Metrics'
## Note: no visible binding for global variable 'caret' 

## set up working directory

#setwd(paste (root, "/myGit/mixturemodel/reconData/para1/", sep=""))
## read in data from txt file
#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))
## read in data from txt file
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")
data <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")
                   
##### BEGIN: data partition >>>>>
## set random seed
set.seed(12345)
#set.seed(34546)
## create data partition

inTrainingSet <- createDataPartition(data$label, p=.7, list=FALSE)
labelTrain <- data[ inTrainingSet,]
labelTest <- data[-inTrainingSet,]
nrow(labelTrain)
nrow(labelTest)


##### END: data partition <<<<<
##### BEGIN: tune the parameters >>>>>
## control:
# resampling technique: 5-repeat 10-fold cross-validation
# performance metrics: ROC AUC curve

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     summaryFunction = multiClassSummary,
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
svmPred <- predict(svmFit, labelTest)
str(svmPred)

## predicted probabilities
svmProbs <- predict(svmFit, labelTest, type = "prob")
str(svmProbs)
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
