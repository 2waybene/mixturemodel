
##  File: predModHong.R
##  Author: Hong Xu
library(caret)
library(pROC)
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

##	param1
#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))
## read in data from txt file

##	param2
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

##	param3
#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")

##	param4
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


setwd(paste (root, "/myGit/mixturemodel/modeling/", sep=""))
sink ("log_param4.txt")

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
#str(svmPred)

## predicted probabilities
svmProbs <- predict(svmFit, labelTest, type = "prob")
#str(svmProbs)
cat ("This is the prediction with SVM")
cat("\n")
confusionMatrix(svmPred, labelTest$label)

##### BEGIN: train model - random forest >>>>>
rfFit <- train(label ~ ., method = "rf", data = labelTrain)
rfPred <- predict(rfFit, labelTest)
cat ("This is the prediction with random forest")
cat("\n")
confusionMatrix(rfPred, labelTest$label)


##### BEGIN: train model - regularized random forest >>>>>
rrfFit <- train(label ~ ., method = "RRF", data = labelTrain)
rrfPred <- predict(rrfFit, labelTest)
cat ("This is the prediction with regularized random forest")
cat("\n")
confusionMatrix(rrfPred, labelTest$label)


##### BEGIN: train model - knn >>>>>
knnFit  <- train(
  label ~ .,
  data = labelTrain,
  method='knn',
  tuneGrid=expand.grid(.k=1:25),
  metric='Accuracy',
  trControl=trainControl(
    method='repeatedcv', 
    number=10, 
    repeats=15))

knnPred <- predict(knnFit , labelTest)
cat ("This is the prediction with k-nearest neighbor")
cat("\n")
confusionMatrix(knnPred, labelTest$label)

##	Neural network

nnetFit <- train(  
	label ~ .,
  	data = labelTrain,
      method = "nnet",
      trace = FALSE,
      maxit = 100)

nnetPred <- predict(nnetFit , labelTest)
cat ("This is the prediction with neural network")
cat("\n")
confusionMatrix(nnetPred, labelTest$label)



##### BEGIN: train model - NaiveBayes >>>>>
nbFit  <- train(
  label ~ .,
  data = labelTrain,
  method='nb',
  trControl=trainControl(method='cv',number=10)
  )

nbPred <- predict(nbFit , labelTest)
cat ("This is the prediction with naive bayes")
cat("\n")
confusionMatrix(nbPred, labelTest$label)

sink()


##### BEGIN: train model - knn3 >>>>>

##	NOT working yet!
knn3Fit  <- knn3(
  label ~ .,
  data = labelTrain,
	k = 11,
  trControl=trainControl(
    method='repeatedcv', 
    number=10, 
    repeats=15))

knn3Pred <- predict(knn3Fit , labelTest)
confusionMatrix(knnPred, labelTest$label)




