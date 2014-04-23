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
## read in data from txt file
#data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))
## read in data from txt file
#data <- read.table("recon_3classes_para2.txt", header=TRUE, sep = "\t")

#data <- read.table("recon_3classes_para3.txt", header=TRUE, sep = "\t")
data <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")


library(caret)
model <- train(
  Species~., 
  data=iris, 
  method='knn',
  tuneGrid=expand.grid(.k=1:25),
  metric='Accuracy',
  trControl=trainControl(
    method='repeatedcv', 
    number=10, 
    repeats=15))

model
plot(model)
confusionMatrix(model)

##========================================================================

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))
oralOncoDX.dt <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")
str(oralOncoDX.dt)

model <- train(
  label~., 
  data=oralOncoDX.dt, 
  method='knn',
  tuneGrid=expand.grid(.k=1:25),
  metric='Accuracy',
  trControl=trainControl(
    method='repeatedcv', 
    number=10, 
    repeats=15))

model
plot(model)
confusionMatrix(model)

